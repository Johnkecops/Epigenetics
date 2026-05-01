"""
tcga_downloader.py
------------------
GDC API wrapper for retrieving TCGA gene expression / miRNA count data
and associated clinical metadata.

Data flow mirrors the TCGAbiolinks GDCquery/GCDdownload/GCDprepare workflow
used in Parikesit et al. (2023) and Fugaha et al. (2024), implemented via
direct REST calls to the NCI Genomic Data Commons API.
"""

import io
import json
import logging
import os
import tarfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import requests

from .constants import (
    CLINICAL_FIELDS,
    DATA_TYPE_CONFIG,
    GDC_BASE_URL,
    GDC_MAX_RETRIES,
    GDC_PAGE_SIZE,
    GDC_RETRY_DELAY,
    NORMAL_CODES,
    TUMOR_CODES,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _gdc_get(endpoint: str, params: dict, retries: int = GDC_MAX_RETRIES) -> dict:
    """GET request against the GDC REST API with retry logic."""
    url = f"{GDC_BASE_URL}/{endpoint}"
    for attempt in range(retries):
        try:
            resp = requests.get(url, params=params, timeout=60)
            resp.raise_for_status()
            return resp.json()
        except requests.RequestException as exc:
            if attempt < retries - 1:
                logger.warning("GDC GET failed (attempt %d): %s. Retrying...", attempt + 1, exc)
                time.sleep(GDC_RETRY_DELAY)
            else:
                raise


def _gdc_post_download(file_ids: List[str]) -> bytes:
    """POST request to GDC /data endpoint to download a batch of files."""
    url = f"{GDC_BASE_URL}/data"
    payload = {"ids": file_ids}
    headers = {"Content-Type": "application/json"}
    resp = requests.post(url, data=json.dumps(payload), headers=headers, timeout=300)
    resp.raise_for_status()
    return resp.content


def _extract_sample_type_code(barcode: str) -> str:
    """
    Extract the two-digit sample type code from a TCGA barcode.
    e.g. TCGA-44-2655-01A-01R-0821-07 -> '01' (primary tumor)
         TCGA-44-2655-11A-01R-0821-07 -> '11' (normal solid tissue)
    """
    parts = barcode.split("-")
    if len(parts) >= 4:
        return parts[3][:2]
    return "99"  # unknown


def _label_sample_type(barcode: str) -> str:
    code = _extract_sample_type_code(barcode)
    if code in TUMOR_CODES:
        return "Tumor"
    if code in NORMAL_CODES:
        return "Normal"
    return "Other"


# ---------------------------------------------------------------------------
# File metadata retrieval
# ---------------------------------------------------------------------------

def query_files(
    project_id: str,
    data_type_key: str,
    progress_callback=None,
) -> pd.DataFrame:
    """
    Query the GDC Files API for all files matching a TCGA project and data type.

    Parameters
    ----------
    project_id : str
        e.g. "TCGA-LUSC"
    data_type_key : str
        Key in DATA_TYPE_CONFIG, e.g. "Gene Expression" or "miRNA Expression"
    progress_callback : callable, optional
        Called with (current, total, message) during pagination.

    Returns
    -------
    pd.DataFrame
        Columns: file_id, file_name, case_id, submitter_id, sample_ids, md5sum
    """
    cfg = DATA_TYPE_CONFIG[data_type_key]

    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": project_id}},
            {"op": "=", "content": {"field": "data_category", "value": cfg["data_category"]}},
            {"op": "=", "content": {"field": "data_type", "value": cfg["data_type"]}},
        ],
    }
    if cfg.get("workflow_type"):
        filters["content"].append(
            {"op": "=", "content": {"field": "analysis.workflow_type", "value": cfg["workflow_type"]}}
        )

    fields = [
        "file_id",
        "file_name",
        "cases.case_id",
        "cases.submitter_id",
        "cases.samples.sample_id",
        "cases.samples.submitter_id",
        "cases.samples.sample_type",
        "cases.samples.sample_type_id",
        "md5sum",
        "file_size",
    ]

    all_hits = []
    from_idx = 0

    # First call to get total count
    params = {
        "filters": json.dumps(filters),
        "fields": ",".join(fields),
        "size": GDC_PAGE_SIZE,
        "from": from_idx,
        "format": "json",
    }
    data = _gdc_get("files", params)
    total = data["data"]["pagination"]["total"]

    if progress_callback:
        progress_callback(0, total, f"Found {total} files. Fetching metadata...")

    all_hits.extend(data["data"]["hits"])
    from_idx += GDC_PAGE_SIZE

    while from_idx < total:
        params["from"] = from_idx
        data = _gdc_get("files", params)
        all_hits.extend(data["data"]["hits"])
        from_idx += GDC_PAGE_SIZE
        if progress_callback:
            progress_callback(min(from_idx, total), total, "Fetching file metadata...")

    # Flatten nested structure
    rows = []
    for hit in all_hits:
        case = hit.get("cases", [{}])[0]
        samples = case.get("samples", [{}])
        for sample in samples:
            rows.append({
                "file_id": hit["file_id"],
                "file_name": hit.get("file_name", ""),
                "file_size": hit.get("file_size", 0),
                "md5sum": hit.get("md5sum", ""),
                "case_id": case.get("case_id", ""),
                "case_submitter_id": case.get("submitter_id", ""),
                "sample_id": sample.get("sample_id", ""),
                "sample_submitter_id": sample.get("submitter_id", ""),
                "sample_type": sample.get("sample_type", ""),
                "sample_type_id": sample.get("sample_type_id", ""),
            })

    return pd.DataFrame(rows).drop_duplicates(subset=["file_id"])


# ---------------------------------------------------------------------------
# Clinical metadata retrieval
# ---------------------------------------------------------------------------

def get_clinical_metadata(
    project_id: str,
    progress_callback=None,
) -> pd.DataFrame:
    """
    Retrieve per-case clinical metadata from the GDC Cases API.

    Returns
    -------
    pd.DataFrame
        Index = case_id; columns = demographic/exposure/diagnosis fields.
    """
    filters = {
        "op": "=",
        "content": {
            "field": "project.project_id",
            "value": project_id,
        },
    }

    all_hits = []
    from_idx = 0

    params = {
        "filters": json.dumps(filters),
        "fields": ",".join(CLINICAL_FIELDS),
        "expand": "demographic,exposures,diagnoses,samples",
        "size": GDC_PAGE_SIZE,
        "from": from_idx,
        "format": "json",
    }
    data = _gdc_get("cases", params)
    total = data["data"]["pagination"]["total"]

    if progress_callback:
        progress_callback(0, total, f"Found {total} cases. Fetching clinical data...")

    all_hits.extend(data["data"]["hits"])
    from_idx += GDC_PAGE_SIZE

    while from_idx < total:
        params["from"] = from_idx
        data = _gdc_get("cases", params)
        all_hits.extend(data["data"]["hits"])
        from_idx += GDC_PAGE_SIZE
        if progress_callback:
            progress_callback(min(from_idx, total), total, "Fetching clinical metadata...")

    rows = []
    for hit in all_hits:
        demo = hit.get("demographic", {}) or {}
        diagnoses = hit.get("diagnoses", [{}])
        diag = diagnoses[0] if diagnoses else {}
        exposures = hit.get("exposures", [{}])
        expo = exposures[0] if exposures else {}

        row = {
            "case_id": hit.get("case_id", ""),
            "case_submitter_id": hit.get("submitter_id", ""),
            "gender": demo.get("gender", "not reported"),
            "race": demo.get("race", "not reported"),
            "ethnicity": demo.get("ethnicity", "not reported"),
            "vital_status": demo.get("vital_status", "not reported"),
            "age_at_index": demo.get("age_at_index"),
            "tobacco_smoking_status": str(expo.get("tobacco_smoking_status", "Unknown")),
            "tobacco_smoking_history": str(expo.get("tobacco_smoking_history_indicator", "Unknown")),
            "pack_years_smoked": expo.get("pack_years_smoked"),
            "primary_diagnosis": diag.get("primary_diagnosis", "not reported"),
            "tumor_stage": diag.get("tumor_stage", diag.get("ajcc_pathologic_stage", "not reported")),
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    df = df.set_index("case_id")
    return df


# ---------------------------------------------------------------------------
# File download and parsing
# ---------------------------------------------------------------------------

def download_and_parse_counts(
    file_df: pd.DataFrame,
    data_type_key: str,
    cache_dir: Optional[str] = None,
    progress_callback=None,
    batch_size: int = 50,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Download count files from GDC and assemble a sample x feature count matrix.

    Files are downloaded in batches and cached locally if cache_dir is provided.

    Parameters
    ----------
    file_df : pd.DataFrame
        Output of query_files(); must have 'file_id', 'sample_submitter_id' columns.
    data_type_key : str
        Key in DATA_TYPE_CONFIG.
    cache_dir : str, optional
        Directory for caching downloaded files. Skips download if cache exists.
    progress_callback : callable, optional
        Called with (current, total, message).
    batch_size : int
        Number of files per download request.

    Returns
    -------
    count_matrix : pd.DataFrame
        Samples x features integer count matrix.
    sample_info : pd.DataFrame
        Sample-level metadata (sample_type, case_id, etc.) aligned to count_matrix.
    """
    cfg = DATA_TYPE_CONFIG[data_type_key]
    cache_path = Path(cache_dir) if cache_dir else None

    file_ids = file_df["file_id"].tolist()
    id_to_barcode = dict(zip(file_df["file_id"], file_df["sample_submitter_id"]))
    id_to_case = dict(zip(file_df["file_id"], file_df["case_id"]))
    id_to_name = dict(zip(file_df["file_id"], file_df["file_name"]))

    total = len(file_ids)
    parsed: Dict[str, pd.Series] = {}

    for batch_start in range(0, total, batch_size):
        batch_ids = file_ids[batch_start: batch_start + batch_size]

        if progress_callback:
            progress_callback(
                batch_start, total,
                f"Downloading files {batch_start + 1}-{min(batch_start + batch_size, total)} of {total}..."
            )

        # Check cache first
        if cache_path:
            cached_all = True
            for fid in batch_ids:
                fname = id_to_name.get(fid, fid)
                if not (cache_path / fname).exists():
                    cached_all = False
                    break
            if cached_all:
                for fid in batch_ids:
                    fname = id_to_name.get(fid, fid)
                    fpath = cache_path / fname
                    series = _parse_count_file(fpath, cfg)
                    barcode = id_to_barcode.get(fid, fid)
                    parsed[barcode] = series
                continue

        # Download batch
        content = _gdc_post_download(batch_ids)
        _extract_and_parse_batch(
            content, batch_ids, id_to_barcode, id_to_name, cfg, parsed, cache_path
        )

    if not parsed:
        raise ValueError("No count data could be parsed from the downloaded files.")

    count_matrix = pd.DataFrame(parsed).T
    count_matrix = count_matrix.fillna(0).astype(int)

    # Build sample_info from file_df aligned to count_matrix index
    sample_info_rows = []
    for barcode in count_matrix.index:
        fid = {v: k for k, v in id_to_barcode.items()}.get(barcode, "")
        case_id = id_to_case.get(fid, "")
        sample_type = _label_sample_type(barcode)
        sample_info_rows.append({
            "sample_id": barcode,
            "case_id": case_id,
            "sample_type": sample_type,
        })

    sample_info = pd.DataFrame(sample_info_rows).set_index("sample_id")

    if progress_callback:
        progress_callback(total, total, "Download complete.")

    return count_matrix, sample_info


def _extract_and_parse_batch(
    content: bytes,
    batch_ids: List[str],
    id_to_barcode: dict,
    id_to_name: dict,
    cfg: dict,
    parsed: dict,
    cache_path: Optional[Path],
):
    """Handle tar archive or single-file response from GDC /data."""
    try:
        # Multi-file response -> tar.gz
        with tarfile.open(fileobj=io.BytesIO(content), mode="r:gz") as tar:
            for member in tar.getmembers():
                if member.isfile():
                    f = tar.extractfile(member)
                    if f is None:
                        continue
                    raw = f.read()
                    # Match to file_id by filename
                    fname = os.path.basename(member.name)
                    fid = _match_filename_to_id(fname, id_to_name)
                    if not fid:
                        continue
                    barcode = id_to_barcode.get(fid, fid)
                    series = _parse_count_bytes(raw, cfg)
                    parsed[barcode] = series
                    if cache_path:
                        cache_path.mkdir(parents=True, exist_ok=True)
                        (cache_path / fname).write_bytes(raw)
    except tarfile.TarError:
        # Single-file response
        if len(batch_ids) == 1:
            fid = batch_ids[0]
            barcode = id_to_barcode.get(fid, fid)
            series = _parse_count_bytes(content, cfg)
            parsed[barcode] = series
            if cache_path:
                cache_path.mkdir(parents=True, exist_ok=True)
                fname = id_to_name.get(fid, fid)
                (cache_path / fname).write_bytes(content)


def _match_filename_to_id(filename: str, id_to_name: dict) -> Optional[str]:
    for fid, fname in id_to_name.items():
        if fname == filename or fname.endswith(filename) or filename.endswith(fname):
            return fid
    return None


def _parse_count_file(path: Path, cfg: dict) -> pd.Series:
    return _parse_count_bytes(path.read_bytes(), cfg)


def _parse_count_bytes(raw: bytes, cfg: dict) -> pd.Series:
    """Parse raw bytes of a count file into a gene->count Series."""
    text = raw.decode("utf-8", errors="replace")
    df = pd.read_csv(io.StringIO(text), sep="\t", comment="#")

    id_col = cfg["id_column"]
    count_col = cfg["count_column"]

    # Drop summary rows present in STAR count files
    skip_prefixes = cfg.get("summary_prefixes", [])
    if skip_prefixes and id_col in df.columns:
        mask = df[id_col].str.startswith(tuple(skip_prefixes))
        df = df[~mask]

    if id_col not in df.columns or count_col not in df.columns:
        # Fallback: assume first col = ID, second = count
        df.columns = [id_col] + list(df.columns[1:])
        if count_col not in df.columns:
            count_col = df.columns[1]

    df = df[[id_col, count_col]].dropna()
    df = df.set_index(id_col)
    return df[count_col].astype(float)


# ---------------------------------------------------------------------------
# Public convenience function
# ---------------------------------------------------------------------------

def fetch_tcga_data(
    project_id: str,
    data_type_key: str,
    cache_dir: Optional[str] = None,
    progress_callback=None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Full data retrieval pipeline: query files -> download counts -> get clinical metadata.

    Returns
    -------
    count_matrix : pd.DataFrame
        Samples x features raw count matrix.
    sample_info : pd.DataFrame
        Sample-level info (sample_type, case_id) aligned to count_matrix rows.
    clinical_df : pd.DataFrame
        Per-case clinical metadata (gender, race, smoking status, etc.).
    """
    if progress_callback:
        progress_callback(0, 4, "Querying GDC file catalog...")

    file_df = query_files(project_id, data_type_key)

    if progress_callback:
        progress_callback(1, 4, f"Found {len(file_df)} files. Fetching clinical metadata...")

    clinical_df = get_clinical_metadata(project_id)

    if progress_callback:
        progress_callback(2, 4, "Downloading count data...")

    count_matrix, sample_info = download_and_parse_counts(
        file_df, data_type_key, cache_dir=cache_dir, progress_callback=progress_callback
    )

    if progress_callback:
        progress_callback(4, 4, "Data retrieval complete.")

    return count_matrix, sample_info, clinical_df
