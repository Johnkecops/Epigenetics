"""
preprocessing.py
----------------
Data preprocessing: clinical filter application, sample-type labeling,
low-expression gene filtering, and DESeq2 input construction.

Implements the preprocessing steps from Parikesit et al. (2023):
- TCGA barcode-based tumor/normal labeling
- Smoking group stratification
- Upper-quartile-inspired expression filtering (>50% samples expressing)
- DESeq2-compatible integer count matrix construction
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .constants import DEFAULT_PARAMS, SMOKING_STATUS_LABELS

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Clinical filter application
# ---------------------------------------------------------------------------

def apply_clinical_filters(
    count_matrix: pd.DataFrame,
    sample_info: pd.DataFrame,
    clinical_df: pd.DataFrame,
    filters: Dict,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Merge clinical data onto sample_info and apply user-selected filters.

    Parameters
    ----------
    count_matrix : pd.DataFrame
        Samples x features raw count matrix.
    sample_info : pd.DataFrame
        Sample-level info with 'case_id' column.
    clinical_df : pd.DataFrame
        Per-case clinical metadata (index = case_id).
    filters : dict
        Keys and selected values, e.g.:
        {
            "gender": ["male"],
            "race": ["white", "black or african american"],
            "tobacco_smoking_status": ["1", "2"],
            "sample_type": ["Tumor", "Normal"],
        }
        Empty list or absent key = no filter for that field.

    Returns
    -------
    filtered_counts : pd.DataFrame
    filtered_metadata : pd.DataFrame
        Combined sample_info + clinical columns, index = sample_id.
    """
    # Merge clinical onto sample_info
    merged = sample_info.copy()
    merged = merged.join(clinical_df, on="case_id", how="left")

    # Apply each filter
    for field, selected in filters.items():
        if not selected:
            continue
        if field not in merged.columns:
            logger.warning("Filter field '%s' not found in metadata. Skipping.", field)
            continue
        selected_lower = [str(s).lower() for s in selected]
        mask = merged[field].astype(str).str.lower().isin(selected_lower)
        merged = merged[mask]

    kept_samples = merged.index.intersection(count_matrix.index)
    filtered_counts = count_matrix.loc[kept_samples]
    filtered_metadata = merged.loc[kept_samples]

    logger.info(
        "After clinical filters: %d samples (from %d).",
        len(filtered_counts), len(count_matrix)
    )
    return filtered_counts, filtered_metadata


# ---------------------------------------------------------------------------
# Expression-level filtering (mirrors paper's >50% sample expression criterion)
# ---------------------------------------------------------------------------

def filter_low_expression(
    count_matrix: pd.DataFrame,
    min_count: int = DEFAULT_PARAMS["min_count"],
    min_samples_fraction: float = DEFAULT_PARAMS["min_samples_fraction"],
) -> pd.DataFrame:
    """
    Remove features with negligible expression.

    Criterion (Parikesit et al. 2023):
      Keep features where at least `min_samples_fraction` of samples have
      raw count >= `min_count`.

    Parameters
    ----------
    count_matrix : pd.DataFrame
        Samples (rows) x features (columns).
    min_count : int
        Minimum count threshold to consider a feature expressed.
    min_samples_fraction : float
        Fraction of samples that must pass min_count.

    Returns
    -------
    pd.DataFrame
        Filtered count matrix.
    """
    n_samples = count_matrix.shape[0]
    min_samples = max(1, int(np.ceil(min_samples_fraction * n_samples)))

    expressed_mask = (count_matrix >= min_count).sum(axis=0) >= min_samples
    filtered = count_matrix.loc[:, expressed_mask]

    logger.info(
        "Low-expression filter: %d/%d features retained (min_count=%d, min_samples=%d).",
        filtered.shape[1], count_matrix.shape[1], min_count, min_samples
    )
    return filtered


# ---------------------------------------------------------------------------
# DESeq2 input construction
# ---------------------------------------------------------------------------

def build_deseq_inputs_tumor_vs_normal(
    count_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: Optional[str] = None,
    group_value: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Construct DESeq2 inputs for tumor vs. normal comparison within an optional subgroup.

    Mirrors DESeq2 Analysis 1 from Parikesit et al. (2023):
      design = ~ tumor
    Used to identify genes aberrantly expressed in tumor development.

    Parameters
    ----------
    count_matrix : pd.DataFrame
        Samples x features.
    metadata : pd.DataFrame
        Must contain 'sample_type' column ('Tumor' or 'Normal').
    group_col : str, optional
        Column to subset on (e.g. 'tobacco_smoking_status').
    group_value : str, optional
        Value to keep within group_col (e.g. '2' for current smokers).

    Returns
    -------
    counts : pd.DataFrame
        Integer counts, samples as rows, features as columns.
    meta : pd.DataFrame
        Metadata with 'condition' column ('Tumor' / 'Normal').
    """
    meta = metadata.copy()

    if group_col and group_value is not None:
        mask = meta[group_col].astype(str) == str(group_value)
        meta = meta[mask]

    meta = meta[meta["sample_type"].isin(["Tumor", "Normal"])]
    counts = count_matrix.loc[meta.index]
    counts = filter_low_expression(counts)
    counts = counts.round().astype(int)

    meta = meta.copy()
    meta["condition"] = meta["sample_type"]

    _validate_groups(meta, "condition", min_samples=3)
    return counts, meta[["condition"]]


def build_deseq_inputs_clinical_variable(
    count_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    clinical_col: str,
    group_a_values: List[str],
    group_b_values: List[str],
    group_a_label: str = "GroupA",
    group_b_label: str = "GroupB",
    include_tumor_covariate: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Construct DESeq2 inputs for comparing two clinical subgroups in tumor samples.

    Mirrors DESeq2 Analysis 2 from Parikesit et al. (2023):
      design = ~ clinical_variable + tumor
    Used to identify genes differentially expressed between (e.g.) smoking groups.

    Parameters
    ----------
    count_matrix : pd.DataFrame
    metadata : pd.DataFrame
    clinical_col : str
        Column to stratify on (e.g. 'tobacco_smoking_status').
    group_a_values : list of str
        Values that define Group A (e.g. ['2', '3', '4'] for smokers).
    group_b_values : list of str
        Values that define Group B (e.g. ['1'] for non-smokers).
    group_a_label : str
    group_b_label : str
    include_tumor_covariate : bool
        Whether to add 'tumor_status' as a blocking factor.

    Returns
    -------
    counts : pd.DataFrame
    meta : pd.DataFrame
        Metadata with 'condition' and optionally 'tumor_status' columns.
    """
    meta = metadata.copy()
    meta[clinical_col] = meta[clinical_col].astype(str)

    mask_a = meta[clinical_col].isin([str(v) for v in group_a_values])
    mask_b = meta[clinical_col].isin([str(v) for v in group_b_values])
    meta = meta[mask_a | mask_b]

    meta["condition"] = np.where(
        meta[clinical_col].isin([str(v) for v in group_a_values]),
        group_a_label,
        group_b_label,
    )

    if include_tumor_covariate and "sample_type" in meta.columns:
        meta["tumor_status"] = np.where(
            meta["sample_type"] == "Tumor", "tumor", "normal"
        )

    counts = count_matrix.loc[meta.index]
    counts = filter_low_expression(counts)
    counts = counts.round().astype(int)

    design_cols = ["condition"]
    if include_tumor_covariate and "tumor_status" in meta.columns:
        design_cols.append("tumor_status")

    _validate_groups(meta, "condition", min_samples=3)
    return counts, meta[design_cols]


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def _validate_groups(meta: pd.DataFrame, condition_col: str, min_samples: int = 3):
    counts = meta[condition_col].value_counts()
    low = counts[counts < min_samples]
    if not low.empty:
        logger.warning(
            "Groups with fewer than %d samples: %s. "
            "DESeq2 may produce unreliable results.",
            min_samples,
            low.to_dict()
        )


def get_available_filters(
    sample_info: pd.DataFrame,
    clinical_df: pd.DataFrame,
) -> Dict[str, List[str]]:
    """
    Return available unique values for each filterable clinical field.

    Used by the Streamlit UI to populate filter widgets dynamically.
    """
    merged = sample_info.join(clinical_df, on="case_id", how="left")

    filterable_cols = [
        "gender", "race", "ethnicity",
        "tobacco_smoking_status", "vital_status",
        "tumor_stage", "sample_type",
    ]

    result = {}
    for col in filterable_cols:
        if col in merged.columns:
            vals = sorted(
                merged[col].dropna().astype(str).unique().tolist()
            )
            # Remove generic 'not reported' from display if other values exist
            if len(vals) > 1 and "not reported" in vals:
                vals = [v for v in vals if v != "not reported"] + ["not reported"]
            result[col] = vals

    return result


def summarize_samples(
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Return a cross-tabulation of sample_type vs. clinical variables for display."""
    cols = [c for c in ["sample_type", "gender", "race", "tobacco_smoking_status"]
            if c in metadata.columns]
    if len(cols) < 2:
        return pd.DataFrame()
    return metadata[cols].value_counts().reset_index(name="count")
