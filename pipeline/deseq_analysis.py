"""
deseq_analysis.py
-----------------
Differential expression analysis using pydeseq2, a Python port of the R DESeq2
package (Love et al., Genome Biology, 2014).

Implements the two-stage analysis pipeline from Parikesit et al. (2023):
  Stage 1 - Tumor vs. Normal within subgroup (design = ~ condition)
  Stage 2 - Clinical variable comparison in tumor/normal context
             (design = ~ condition + tumor_status)
  Stage 3 - Overlap DEGs from both stages to find dual-susceptible genes
"""

import logging
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _import_pydeseq2():
    """Lazy import to give a clear error if pydeseq2 is not installed."""
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        return DeseqDataSet, DeseqStats
    except ImportError as exc:
        raise ImportError(
            "pydeseq2 is required for differential expression analysis. "
            "Install with: pip install pydeseq2"
        ) from exc


# ---------------------------------------------------------------------------
# Core DESeq2 wrapper
# ---------------------------------------------------------------------------

def run_deseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    contrast_col: str,
    numerator: str,
    denominator: str,
    design_factors: Optional[List[str]] = None,
    lfc_threshold: float = 0.0,
    quiet: bool = True,
) -> pd.DataFrame:
    """
    Run DESeq2 differential expression analysis.

    Parameters
    ----------
    counts : pd.DataFrame
        Samples (rows) x features (columns), integer counts.
    metadata : pd.DataFrame
        Must have index matching counts.index; must contain contrast_col.
    contrast_col : str
        Column in metadata to contrast (e.g. 'condition').
    numerator : str
        The 'treatment' level (e.g. 'Tumor' or 'GroupA').
    denominator : str
        The 'reference' level (e.g. 'Normal' or 'GroupB').
    design_factors : list of str, optional
        All factors in the design formula. If None, uses [contrast_col].
    lfc_threshold : float
        Log2 fold change threshold for lfcShrink (default 0 = report all).
    quiet : bool
        Suppress pydeseq2 verbose output.

    Returns
    -------
    pd.DataFrame
        Results table with columns:
        baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
    """
    DeseqDataSet, DeseqStats = _import_pydeseq2()

    if design_factors is None:
        design_factors = [contrast_col]

    # Ensure reference level is first (DESeq2 uses alphabetical or explicit)
    metadata = metadata.copy()
    metadata[contrast_col] = pd.Categorical(
        metadata[contrast_col],
        categories=[denominator, numerator],
        ordered=False,
    )

    # Remove samples not in either group
    valid_mask = metadata[contrast_col].isin([numerator, denominator])
    metadata = metadata[valid_mask]
    counts = counts.loc[metadata.index]

    # Ensure integer counts
    counts = counts.astype(int)

    with warnings.catch_warnings():
        if quiet:
            warnings.simplefilter("ignore")

        dds = DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design_factors=design_factors,
            refit_cooks=True,
            quiet=quiet,
        )
        dds.deseq2()

        stat_res = DeseqStats(
            dds,
            contrast=[contrast_col, numerator, denominator],
            quiet=quiet,
        )
        stat_res.summary()

    results = stat_res.results_df.copy()
    results.index.name = "feature_id"
    results = results.sort_values("padj")
    return results


# ---------------------------------------------------------------------------
# Filtering helpers
# ---------------------------------------------------------------------------

def filter_degs(
    results: pd.DataFrame,
    lfc_threshold: float,
    padj_threshold: float,
) -> pd.DataFrame:
    """
    Apply LFC and adjusted p-value thresholds to extract significant DEGs.

    Criterion from Parikesit et al. (2023):
    - |log2FoldChange| >= lfc_threshold
    - padj < padj_threshold
    """
    mask = (
        (results["padj"] < padj_threshold) &
        (results["log2FoldChange"].abs() >= lfc_threshold)
    )
    return results[mask].copy()


def annotate_expression_direction(results: pd.DataFrame) -> pd.DataFrame:
    """Add 'expression' column: 'Overexpressed', 'Underexpressed', 'Unchanged'."""
    results = results.copy()
    results["expression"] = "Unchanged"
    results.loc[results["log2FoldChange"] > 0, "expression"] = "Overexpressed"
    results.loc[results["log2FoldChange"] < 0, "expression"] = "Underexpressed"
    return results


# ---------------------------------------------------------------------------
# Two-stage pipeline with gene overlap
# ---------------------------------------------------------------------------

def run_dual_stage_pipeline(
    counts_stage1: pd.DataFrame,
    meta_stage1: pd.DataFrame,
    counts_stage2: pd.DataFrame,
    meta_stage2: pd.DataFrame,
    stage1_numerator: str = "Tumor",
    stage1_denominator: str = "Normal",
    stage2_numerator: str = "GroupA",
    stage2_denominator: str = "GroupB",
    lfc_threshold: float = 1.0,
    padj_threshold: float = 0.01,
    overlap_lfc_threshold: float = 2.0,
    progress_callback=None,
) -> Dict:
    """
    Full two-stage differential expression pipeline mirroring Parikesit et al. (2023).

    Stage 1: Tumor vs. Normal within clinical subgroup
             -> identifies tumor-associated DEGs
    Stage 2: Clinical variable comparison (e.g. smoker vs. non-smoker)
             -> identifies clinically susceptible DEGs
    Overlap: Intersection of Stage 1 and Stage 2 DEG lists
             -> dual-susceptible biomarker candidates

    Returns
    -------
    dict with keys:
        results_stage1 : full DESeq2 results table (Stage 1)
        results_stage2 : full DESeq2 results table (Stage 2)
        degs_stage1    : significant DEGs from Stage 1
        degs_stage2    : significant DEGs from Stage 2
        overlap        : DEGs in both stages, filtered by overlap_lfc_threshold
        summary        : dict of counts for reporting
    """
    if progress_callback:
        progress_callback(0, 4, "Running Stage 1: Tumor vs. Normal DESeq2 analysis...")

    results_stage1 = run_deseq2(
        counts_stage1, meta_stage1,
        contrast_col="condition",
        numerator=stage1_numerator,
        denominator=stage1_denominator,
        design_factors=["condition"],
    )

    if progress_callback:
        progress_callback(1, 4, "Running Stage 2: Clinical group DESeq2 analysis...")

    design_factors_2 = list(meta_stage2.columns)
    results_stage2 = run_deseq2(
        counts_stage2, meta_stage2,
        contrast_col="condition",
        numerator=stage2_numerator,
        denominator=stage2_denominator,
        design_factors=design_factors_2,
    )

    if progress_callback:
        progress_callback(2, 4, "Filtering DEGs and computing overlap...")

    degs_stage1 = filter_degs(results_stage1, lfc_threshold, padj_threshold)
    degs_stage2 = filter_degs(results_stage2, lfc_threshold, padj_threshold)

    degs_stage1 = annotate_expression_direction(degs_stage1)
    degs_stage2 = annotate_expression_direction(degs_stage2)

    overlap_genes = set(degs_stage1.index) & set(degs_stage2.index)

    # Build overlap table filtered by stricter LFC threshold (paper uses |LFC| >= 2)
    if overlap_genes:
        overlap_df = results_stage1.loc[
            results_stage1.index.isin(overlap_genes)
        ].copy()
        overlap_df = overlap_df[
            overlap_df["log2FoldChange"].abs() >= overlap_lfc_threshold
        ]
        overlap_df = overlap_df.sort_values("log2FoldChange")
        overlap_df = annotate_expression_direction(overlap_df)
    else:
        overlap_df = pd.DataFrame()

    if progress_callback:
        progress_callback(4, 4, "Analysis complete.")

    summary = {
        "n_stage1_tested": len(results_stage1),
        "n_stage1_degs": len(degs_stage1),
        "n_stage2_tested": len(results_stage2),
        "n_stage2_degs": len(degs_stage2),
        "n_overlap_raw": len(overlap_genes),
        "n_overlap_filtered": len(overlap_df) if not overlap_df.empty else 0,
        "stage1_overexpressed": int((degs_stage1["log2FoldChange"] > 0).sum()),
        "stage1_underexpressed": int((degs_stage1["log2FoldChange"] < 0).sum()),
        "stage2_overexpressed": int((degs_stage2["log2FoldChange"] > 0).sum()),
        "stage2_underexpressed": int((degs_stage2["log2FoldChange"] < 0).sum()),
    }

    return {
        "results_stage1": results_stage1,
        "results_stage2": results_stage2,
        "degs_stage1": degs_stage1,
        "degs_stage2": degs_stage2,
        "overlap": overlap_df,
        "summary": summary,
    }


# ---------------------------------------------------------------------------
# Single-stage pipeline (e.g. miRNA study - Fugaha et al. 2024)
# ---------------------------------------------------------------------------

def run_single_stage_pipeline(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    numerator: str = "Tumor",
    denominator: str = "Normal",
    lfc_threshold: float = 1.0,
    padj_threshold: float = 0.01,
    progress_callback=None,
) -> Dict:
    """
    Single-stage DESeq2 pipeline for tumor vs. normal comparison.

    Mirrors the pipeline in Fugaha et al. (2024) for miRNA expression in
    colorectal carcinoma.

    Returns
    -------
    dict with keys:
        results   : full DESeq2 results table
        degs      : significant DEGs after threshold filtering
        summary   : dict of counts
    """
    if progress_callback:
        progress_callback(0, 2, "Running DESeq2 analysis...")

    results = run_deseq2(
        counts, metadata,
        contrast_col="condition",
        numerator=numerator,
        denominator=denominator,
        design_factors=["condition"],
    )

    if progress_callback:
        progress_callback(1, 2, "Filtering significant features...")

    degs = filter_degs(results, lfc_threshold, padj_threshold)
    degs = annotate_expression_direction(degs)
    results = annotate_expression_direction(results)

    summary = {
        "n_tested": len(results),
        "n_degs": len(degs),
        "n_overexpressed": int((degs["log2FoldChange"] > 0).sum()),
        "n_underexpressed": int((degs["log2FoldChange"] < 0).sum()),
    }

    if progress_callback:
        progress_callback(2, 2, "Analysis complete.")

    return {
        "results": results,
        "degs": degs,
        "summary": summary,
    }


# ---------------------------------------------------------------------------
# Export helpers
# ---------------------------------------------------------------------------

def get_top_features(
    results: pd.DataFrame,
    n: int = 10,
    sort_by: str = "padj",
    direction: Optional[str] = None,
) -> pd.DataFrame:
    """
    Extract top N features from a DESeq2 results table.

    Parameters
    ----------
    direction : str or None
        'up', 'down', or None (both).
    """
    df = results.dropna(subset=["padj", "log2FoldChange"]).copy()

    if direction == "up":
        df = df[df["log2FoldChange"] > 0]
    elif direction == "down":
        df = df[df["log2FoldChange"] < 0]

    return df.sort_values(sort_by, ascending=(sort_by == "padj")).head(n)
