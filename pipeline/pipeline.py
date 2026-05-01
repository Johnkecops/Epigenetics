"""
pipeline.py
-----------
Orchestrates the full TCGA cancer epigenomics analysis pipeline.

Combines data retrieval, preprocessing, DESeq2 analysis, and
visualization into a single callable interface.

Based on methodology from:
  Parikesit et al. (Berita Biologi, 2023) - Gene expression in LUSC
  Fugaha et al. (ICIC, 2024) - miRNA expression in colorectal carcinoma
"""

import logging
import os
from pathlib import Path
from typing import Callable, Dict, List, Optional

import pandas as pd

from .constants import DEFAULT_PARAMS, SMOKING_STATUS_LABELS
from .deseq_analysis import (
    run_dual_stage_pipeline,
    run_single_stage_pipeline,
)
from .preprocessing import (
    apply_clinical_filters,
    build_deseq_inputs_clinical_variable,
    build_deseq_inputs_tumor_vs_normal,
    get_available_filters,
    summarize_samples,
)
from .tcga_downloader import fetch_tcga_data

logger = logging.getLogger(__name__)


class EpigenomicsPipeline:
    """
    End-to-end TCGA cancer epigenomics differential expression pipeline.

    Usage
    -----
    pipeline = EpigenomicsPipeline(cache_dir="./cache")
    pipeline.load_data("TCGA-LUSC", "Gene Expression")
    pipeline.set_clinical_filters({"gender": ["male"]})
    results = pipeline.run_analysis(
        analysis_mode="dual_stage",
        stage2_col="tobacco_smoking_status",
        stage2_group_a=["2", "3", "4"],
        stage2_group_b=["1"],
    )
    """

    def __init__(self, cache_dir: Optional[str] = None):
        self.cache_dir = cache_dir
        self.count_matrix: Optional[pd.DataFrame] = None
        self.sample_info: Optional[pd.DataFrame] = None
        self.clinical_df: Optional[pd.DataFrame] = None
        self.filtered_counts: Optional[pd.DataFrame] = None
        self.filtered_metadata: Optional[pd.DataFrame] = None
        self.results: Optional[Dict] = None
        self._project_id: Optional[str] = None
        self._data_type_key: Optional[str] = None

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def load_data(
        self,
        project_id: str,
        data_type_key: str,
        progress_callback: Optional[Callable] = None,
    ) -> None:
        """
        Download (or load from cache) TCGA count data and clinical metadata.

        Parameters
        ----------
        project_id : str
            e.g. "TCGA-LUSC"
        data_type_key : str
            "Gene Expression" or "miRNA Expression"
        progress_callback : callable, optional
            (current, total, message) -> None
        """
        self._project_id = project_id
        self._data_type_key = data_type_key

        cache_subdir = None
        if self.cache_dir:
            cache_subdir = os.path.join(
                self.cache_dir,
                project_id,
                data_type_key.replace(" ", "_"),
            )
            Path(cache_subdir).mkdir(parents=True, exist_ok=True)

        self.count_matrix, self.sample_info, self.clinical_df = fetch_tcga_data(
            project_id=project_id,
            data_type_key=data_type_key,
            cache_dir=cache_subdir,
            progress_callback=progress_callback,
        )

        logger.info(
            "Loaded %d samples x %d features for %s (%s).",
            self.count_matrix.shape[0],
            self.count_matrix.shape[1],
            project_id,
            data_type_key,
        )

    @property
    def available_filters(self) -> Dict[str, List[str]]:
        """Available values for each filterable clinical field."""
        if self.sample_info is None or self.clinical_df is None:
            return {}
        return get_available_filters(self.sample_info, self.clinical_df)

    @property
    def sample_summary(self) -> pd.DataFrame:
        """Cross-tabulation of sample types and key clinical variables."""
        if self.filtered_metadata is not None:
            return summarize_samples(self.filtered_metadata)
        if self.sample_info is not None and self.clinical_df is not None:
            merged = self.sample_info.join(self.clinical_df, on="case_id", how="left")
            return summarize_samples(merged)
        return pd.DataFrame()

    # ------------------------------------------------------------------
    # Clinical filtering
    # ------------------------------------------------------------------

    def set_clinical_filters(self, filters: Dict) -> None:
        """
        Apply clinical filters to restrict analysis to a patient subgroup.

        Parameters
        ----------
        filters : dict
            e.g. {"gender": ["male"], "tobacco_smoking_status": ["2", "3", "4"]}
            Empty list or missing key = include all.
        """
        if self.count_matrix is None:
            raise RuntimeError("Call load_data() before set_clinical_filters().")

        self.filtered_counts, self.filtered_metadata = apply_clinical_filters(
            self.count_matrix,
            self.sample_info,
            self.clinical_df,
            filters,
        )

        logger.info(
            "Clinical filters applied: %d samples retained.",
            len(self.filtered_counts)
        )

    # ------------------------------------------------------------------
    # Analysis
    # ------------------------------------------------------------------

    def run_analysis(
        self,
        analysis_mode: str = "single_stage",
        stage2_col: Optional[str] = None,
        stage2_group_a: Optional[List[str]] = None,
        stage2_group_b: Optional[List[str]] = None,
        stage2_label_a: str = "Group A",
        stage2_label_b: str = "Group B",
        lfc_threshold: float = DEFAULT_PARAMS["lfc_threshold"],
        padj_threshold: float = DEFAULT_PARAMS["padj_threshold"],
        overlap_lfc_threshold: float = DEFAULT_PARAMS["overlap_lfc_threshold"],
        progress_callback: Optional[Callable] = None,
    ) -> Dict:
        """
        Run differential expression analysis.

        Parameters
        ----------
        analysis_mode : str
            "single_stage" - tumor vs. normal only (Fugaha et al. 2024 style).
            "dual_stage"   - tumor vs. normal + clinical variable comparison
                             with gene overlap (Parikesit et al. 2023 style).
        stage2_col : str, optional
            Clinical column for Stage 2 comparison (required for dual_stage).
        stage2_group_a : list of str
            Values defining Group A in stage2_col.
        stage2_group_b : list of str
            Values defining Group B in stage2_col.
        stage2_label_a : str
        stage2_label_b : str
        lfc_threshold : float
        padj_threshold : float
        overlap_lfc_threshold : float
            Stricter LFC for overlap table (paper uses |LFC| >= 2).
        progress_callback : callable, optional

        Returns
        -------
        dict
            Full analysis results including DESeq2 tables, DEG lists,
            and (for dual_stage) overlap gene table.
        """
        counts = self.filtered_counts if self.filtered_counts is not None else self.count_matrix
        metadata = self.filtered_metadata if self.filtered_metadata is not None else (
            self.sample_info.join(self.clinical_df, on="case_id", how="left")
            if self.sample_info is not None and self.clinical_df is not None
            else self.sample_info
        )

        if counts is None or metadata is None:
            raise RuntimeError("Call load_data() before run_analysis().")

        if analysis_mode == "single_stage":
            counts_in, meta_in = build_deseq_inputs_tumor_vs_normal(counts, metadata)
            self.results = run_single_stage_pipeline(
                counts_in, meta_in,
                numerator="Tumor",
                denominator="Normal",
                lfc_threshold=lfc_threshold,
                padj_threshold=padj_threshold,
                progress_callback=progress_callback,
            )

        elif analysis_mode == "dual_stage":
            if not stage2_col or not stage2_group_a or not stage2_group_b:
                raise ValueError(
                    "dual_stage mode requires stage2_col, stage2_group_a, stage2_group_b."
                )

            # Stage 1: Tumor vs Normal in Group A patients
            counts_s1, meta_s1 = build_deseq_inputs_tumor_vs_normal(
                counts, metadata,
                group_col=stage2_col,
                group_value=stage2_group_a[0] if len(stage2_group_a) == 1 else None,
            )

            # Stage 2: Group A vs Group B (accounting for tumor status)
            counts_s2, meta_s2 = build_deseq_inputs_clinical_variable(
                counts, metadata,
                clinical_col=stage2_col,
                group_a_values=stage2_group_a,
                group_b_values=stage2_group_b,
                group_a_label=stage2_label_a,
                group_b_label=stage2_label_b,
                include_tumor_covariate=True,
            )

            self.results = run_dual_stage_pipeline(
                counts_s1, meta_s1,
                counts_s2, meta_s2,
                stage1_numerator="Tumor",
                stage1_denominator="Normal",
                stage2_numerator=stage2_label_a,
                stage2_denominator=stage2_label_b,
                lfc_threshold=lfc_threshold,
                padj_threshold=padj_threshold,
                overlap_lfc_threshold=overlap_lfc_threshold,
                progress_callback=progress_callback,
            )

        else:
            raise ValueError(f"Unknown analysis_mode: '{analysis_mode}'. "
                             "Choose 'single_stage' or 'dual_stage'.")

        return self.results

    # ------------------------------------------------------------------
    # Convenience properties
    # ------------------------------------------------------------------

    @property
    def deg_table(self) -> Optional[pd.DataFrame]:
        """Main DEG results table from the most recent analysis."""
        if self.results is None:
            return None
        return self.results.get("degs") or self.results.get("degs_stage1")

    @property
    def overlap_table(self) -> Optional[pd.DataFrame]:
        """Overlap DEG table (dual_stage mode only)."""
        if self.results is None:
            return None
        return self.results.get("overlap")

    @property
    def analysis_summary(self) -> Optional[Dict]:
        """Summary statistics from the most recent analysis."""
        if self.results is None:
            return None
        return self.results.get("summary")
