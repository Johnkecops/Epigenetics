"""
app.py
------
Streamlit front-end for the TCGA Cancer Epigenomics Pipeline.

Implements an interactive, step-by-step workflow for:
1. Selecting cancer type and molecular data layer
2. Applying clinical filters (gender, ethnicity, smoking status, etc.)
3. Configuring analysis parameters and comparison groups
4. Running DESeq2 differential expression analysis
5. Visualizing and exporting results

Author: Dr. Arli Aditya Parikesit (arli.parikesit@i3l.ac.id)
Institution: Department of Bioinformatics, i3L University, Jakarta, Indonesia
"""

import logging
import os

import pandas as pd
import streamlit as st

from pipeline.constants import (
    CANCER_TYPES,
    DATA_TYPE_CONFIG,
    DEFAULT_PARAMS,
    SMOKING_STATUS_LABELS,
)
from pipeline.pipeline import EpigenomicsPipeline
from pipeline.visualization import (
    plot_expression_heatmap,
    plot_ma,
    plot_sample_distribution,
    plot_top_features,
    plot_venn,
    plot_volcano,
    plot_wald_scatter,
)

logging.basicConfig(level=logging.INFO)

# ---------------------------------------------------------------------------
# Page configuration
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="TCGA Cancer Epigenomics Pipeline",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("TCGA Cancer Epigenomics Pipeline")
st.markdown(
    """
    Differential gene and miRNA expression analysis using TCGA data.
    Based on methodology from [Parikesit et al. (Berita Biologi, 2023)](https://doi.org/10.14203/beritabiology.v20i1.3991)
    and [Fugaha et al. (ICIC, 2024)](https://doi.org/10.1109/ICIC62776.2024).

    **Workflow:** Select cancer → Filter patients → Configure analysis → Run → Explore results.
    """,
    unsafe_allow_html=False,
)

# ---------------------------------------------------------------------------
# Session state initialization
# ---------------------------------------------------------------------------

def _init_state():
    defaults = {
        "pipeline": None,
        "data_loaded": False,
        "filters_applied": False,
        "analysis_done": False,
        "available_filters": {},
        "results": None,
        "project_id": None,
        "data_type_key": None,
    }
    for key, val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = val

_init_state()

# ---------------------------------------------------------------------------
# Sidebar: Step 1 - Cancer and data type selection
# ---------------------------------------------------------------------------

with st.sidebar:
    st.header("Step 1: Data Selection")

    cancer_labels = {
        f"{code} - {name}": code
        for code, name in CANCER_TYPES.items()
    }
    selected_label = st.selectbox(
        "Cancer Type",
        list(cancer_labels.keys()),
        index=list(cancer_labels.keys()).index("TCGA-LUSC - Lung Squamous Cell Carcinoma"),
        help="Select a TCGA cancer project to analyze.",
    )
    project_id = cancer_labels[selected_label]

    data_type_key = st.selectbox(
        "Data Layer",
        list(DATA_TYPE_CONFIG.keys()),
        help="Gene Expression (mRNA counts) or miRNA Expression.",
    )

    cache_dir = st.text_input(
        "Local Cache Directory",
        value="./tcga_cache",
        help="Downloaded files are cached here to avoid re-downloading.",
    )

    load_btn = st.button("Load Data from TCGA", type="primary", use_container_width=True)

    if load_btn:
        st.session_state.pipeline = EpigenomicsPipeline(cache_dir=cache_dir or None)
        st.session_state.data_loaded = False
        st.session_state.filters_applied = False
        st.session_state.analysis_done = False
        st.session_state.project_id = project_id
        st.session_state.data_type_key = data_type_key

        progress_bar = st.progress(0, text="Initializing...")

        def _progress(current, total, msg):
            pct = int(current / max(total, 1) * 100)
            progress_bar.progress(pct, text=msg)

        try:
            st.session_state.pipeline.load_data(
                project_id, data_type_key, progress_callback=_progress
            )
            st.session_state.available_filters = (
                st.session_state.pipeline.available_filters
            )
            st.session_state.data_loaded = True
            progress_bar.empty()
            st.success(
                f"Loaded {st.session_state.pipeline.count_matrix.shape[0]} samples "
                f"x {st.session_state.pipeline.count_matrix.shape[1]} features."
            )
        except Exception as exc:
            progress_bar.empty()
            st.error(f"Data loading failed: {exc}")

    st.divider()

# ---------------------------------------------------------------------------
# Sidebar: Step 2 - Clinical Filters
# ---------------------------------------------------------------------------

    st.header("Step 2: Clinical Filters")
    st.caption("Leave empty to include all patients.")

    filters = {}
    avail = st.session_state.available_filters

    FILTER_LABELS = {
        "sample_type": "Sample Type",
        "gender": "Gender",
        "race": "Race / Ethnicity",
        "ethnicity": "Hispanic/Latino Ethnicity",
        "tobacco_smoking_status": "Smoking Status",
        "vital_status": "Vital Status",
        "tumor_stage": "Tumor Stage",
    }

    SMOKING_HELP = (
        "TCGA codes: 1=Non-smoker, 2=Current, 3=Reformed >15yr, "
        "4=Reformed <=15yr, 5=Reformed (unknown duration)"
    )

    for col, label in FILTER_LABELS.items():
        if col in avail and avail[col]:
            opts = avail[col]
            help_text = SMOKING_HELP if col == "tobacco_smoking_status" else None
            selected = st.multiselect(label, opts, default=[], help=help_text)
            if selected:
                filters[col] = selected

    apply_filters_btn = st.button(
        "Apply Filters",
        disabled=not st.session_state.data_loaded,
        use_container_width=True,
    )

    if apply_filters_btn:
        try:
            st.session_state.pipeline.set_clinical_filters(filters)
            st.session_state.filters_applied = True
            n = len(st.session_state.pipeline.filtered_counts)
            st.success(f"{n} samples after filtering.")
        except Exception as exc:
            st.error(f"Filter error: {exc}")

    st.divider()

# ---------------------------------------------------------------------------
# Sidebar: Step 3 - Analysis Configuration
# ---------------------------------------------------------------------------

    st.header("Step 3: Analysis Parameters")

    analysis_mode = st.radio(
        "Analysis Mode",
        ["single_stage", "dual_stage"],
        format_func=lambda x: {
            "single_stage": "Tumor vs. Normal (single comparison)",
            "dual_stage": "Dual comparison + gene overlap",
        }[x],
        help=(
            "Single: tumor vs. normal only (Fugaha et al. 2024 style).\n"
            "Dual: tumor/normal + clinical variable with gene overlap "
            "(Parikesit et al. 2023 style)."
        ),
    )

    stage2_col = stage2_group_a = stage2_group_b = None
    stage2_label_a = stage2_label_b = ""

    if analysis_mode == "dual_stage":
        st.subheader("Stage 2 Clinical Comparison")
        clinical_cols = [
            c for c in ["tobacco_smoking_status", "gender", "race",
                        "ethnicity", "tumor_stage", "vital_status"]
            if c in avail
        ]
        stage2_col = st.selectbox(
            "Stratification Variable",
            clinical_cols if clinical_cols else ["tobacco_smoking_status"],
            help="Clinical variable to compare between two groups.",
        )

        if stage2_col and stage2_col in avail:
            col_opts = avail[stage2_col]
            stage2_group_a = st.multiselect(
                "Group A (case/condition)",
                col_opts,
                default=col_opts[:1] if col_opts else [],
                help="e.g. Current smokers: select '2'",
            )
            stage2_group_b = st.multiselect(
                "Group B (reference)",
                col_opts,
                default=col_opts[-1:] if col_opts else [],
                help="e.g. Non-smokers: select '1'",
            )
            stage2_label_a = st.text_input("Group A Label", value="Group A")
            stage2_label_b = st.text_input("Group B Label", value="Group B")

    st.subheader("DESeq2 Thresholds")
    lfc_threshold = st.slider(
        "Log2 Fold Change threshold",
        min_value=0.5, max_value=5.0, value=1.0, step=0.25,
        help="Minimum |LFC| to call a DEG.",
    )
    padj_threshold = st.select_slider(
        "Adjusted P-value threshold",
        options=[0.001, 0.005, 0.01, 0.05, 0.1],
        value=0.01,
        help="FDR-adjusted p-value cutoff.",
    )
    overlap_lfc_threshold = st.slider(
        "Overlap LFC threshold",
        min_value=1.0, max_value=5.0, value=2.0, step=0.5,
        help="Stricter LFC for the final overlap/biomarker table.",
    )

    run_btn = st.button(
        "Run Analysis",
        type="primary",
        disabled=not st.session_state.data_loaded,
        use_container_width=True,
    )

# ---------------------------------------------------------------------------
# Main panel: orchestrate analysis run
# ---------------------------------------------------------------------------

if run_btn:
    if not st.session_state.filters_applied:
        st.session_state.pipeline.set_clinical_filters({})
        st.session_state.filters_applied = True

    progress_placeholder = st.empty()
    prog_bar = progress_placeholder.progress(0, text="Starting analysis...")

    def _run_progress(current, total, msg):
        pct = int(current / max(total, 1) * 100)
        prog_bar.progress(pct, text=msg)

    try:
        results = st.session_state.pipeline.run_analysis(
            analysis_mode=analysis_mode,
            stage2_col=stage2_col,
            stage2_group_a=stage2_group_a if stage2_group_a else None,
            stage2_group_b=stage2_group_b if stage2_group_b else None,
            stage2_label_a=stage2_label_a or "Group A",
            stage2_label_b=stage2_label_b or "Group B",
            lfc_threshold=lfc_threshold,
            padj_threshold=padj_threshold,
            overlap_lfc_threshold=overlap_lfc_threshold,
            progress_callback=_run_progress,
        )
        st.session_state.results = results
        st.session_state.analysis_done = True
        progress_placeholder.empty()
        st.success("Analysis complete.")
    except Exception as exc:
        progress_placeholder.empty()
        st.error(f"Analysis failed: {exc}")
        st.exception(exc)

# ---------------------------------------------------------------------------
# Results display
# ---------------------------------------------------------------------------

if not st.session_state.data_loaded:
    st.info("Configure settings in the sidebar and click **Load Data from TCGA** to begin.")

elif not st.session_state.analysis_done:
    # Show data overview before analysis
    st.subheader("Data Overview")
    pl = st.session_state.pipeline
    col1, col2, col3 = st.columns(3)
    col1.metric("Total Samples", pl.count_matrix.shape[0] if pl.count_matrix is not None else 0)
    col2.metric("Total Features", pl.count_matrix.shape[1] if pl.count_matrix is not None else 0)
    col3.metric(
        "After Filters",
        len(pl.filtered_counts) if pl.filtered_counts is not None else pl.count_matrix.shape[0]
    )

    if pl.sample_info is not None and pl.clinical_df is not None:
        meta = pl.sample_info.join(pl.clinical_df, on="case_id", how="left")
        if "sample_type" in meta.columns:
            st.plotly_chart(
                plot_sample_distribution(meta, "sample_type", "Sample Type Distribution"),
                use_container_width=True
            )
        if "gender" in meta.columns:
            st.plotly_chart(
                plot_sample_distribution(meta, "gender", "Gender Distribution"),
                use_container_width=True
            )

    st.info("Configure analysis parameters in the sidebar and click **Run Analysis**.")

else:
    results = st.session_state.results
    pl = st.session_state.pipeline
    is_dual = analysis_mode == "dual_stage"

    # ----------------------------------------------------------------
    # Summary metrics
    # ----------------------------------------------------------------
    st.subheader("Analysis Summary")
    summary = results.get("summary", {})

    if is_dual:
        col1, col2, col3, col4, col5 = st.columns(5)
        col1.metric("Stage 1 DEGs", summary.get("n_stage1_degs", 0))
        col2.metric("Stage 2 DEGs", summary.get("n_stage2_degs", 0))
        col3.metric("Overlap (raw)", summary.get("n_overlap_raw", 0))
        col4.metric(f"Overlap |LFC|≥{overlap_lfc_threshold}", summary.get("n_overlap_filtered", 0))
        col5.metric("Stage 1 Tested", summary.get("n_stage1_tested", 0))
    else:
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Tested Features", summary.get("n_tested", 0))
        col2.metric("Significant DEGs", summary.get("n_degs", 0))
        col3.metric("Overexpressed", summary.get("n_overexpressed", 0))
        col4.metric("Underexpressed", summary.get("n_underexpressed", 0))

    # ----------------------------------------------------------------
    # Tabs
    # ----------------------------------------------------------------
    tab_names = ["Volcano Plot", "MA Plot", "Top Features", "DEG Table"]
    if is_dual:
        tab_names += ["Venn Diagram", "Wald Scatter", "Overlap Table"]
    tab_names += ["Sample Overview"]

    tabs = st.tabs(tab_names)
    tab_idx = 0

    # -- Volcano --
    with tabs[tab_idx]:
        tab_idx += 1
        st.subheader("Volcano Plot")

        if is_dual:
            stage_choice = st.radio(
                "Select stage", ["Stage 1 (Tumor vs. Normal)", "Stage 2 (Clinical comparison)"],
                horizontal=True, key="volcano_stage"
            )
            r_key = "results_stage1" if "Stage 1" in stage_choice else "results_stage2"
            deg_key = "degs_stage1" if "Stage 1" in stage_choice else "degs_stage2"
            res = results[r_key]
            title = f"Volcano Plot - {stage_choice.split('(')[1].rstrip(')')}"
        else:
            res = results["results"]
            title = "Volcano Plot"
            deg_key = "degs"

        top_genes = list(results[deg_key].head(5).index) if not results[deg_key].empty else []
        fig = plot_volcano(res, lfc_threshold, padj_threshold, highlight_genes=top_genes, title=title)
        st.plotly_chart(fig, use_container_width=True)
        st.caption(
            f"Blue points: significant DEGs (|LFC| >= {lfc_threshold}, padj < {padj_threshold}). "
            "Vertical red lines: LFC thresholds. Horizontal red line: padj threshold."
        )

    # -- MA Plot --
    with tabs[tab_idx]:
        tab_idx += 1
        st.subheader("MA Plot")

        if is_dual:
            stage_choice = st.radio(
                "Select stage", ["Stage 1 (Tumor vs. Normal)", "Stage 2 (Clinical comparison)"],
                horizontal=True, key="ma_stage"
            )
            r_key = "results_stage1" if "Stage 1" in stage_choice else "results_stage2"
            deg_key = "degs_stage1" if "Stage 1" in stage_choice else "degs_stage2"
            res = results[r_key]
            title = f"MA Plot - {stage_choice.split('(')[1].rstrip(')')}"
        else:
            res = results["results"]
            title = "MA Plot"
            deg_key = "degs"

        top_genes = list(results[deg_key].head(5).index) if not results[deg_key].empty else []
        fig = plot_ma(res, lfc_threshold, padj_threshold, highlight_genes=top_genes, title=title)
        st.plotly_chart(fig, use_container_width=True)
        st.caption(
            "X-axis: log10(mean normalized counts). Y-axis: log2 fold change. "
            "Blue points: significant DEGs."
        )

    # -- Top Features --
    with tabs[tab_idx]:
        tab_idx += 1
        st.subheader("Top Differentially Expressed Features")

        deg_df = results.get("degs") or results.get("degs_stage1", pd.DataFrame())
        top_n = st.slider("Number of features to display", 5, 50, 20, key="top_n")

        if not deg_df.empty:
            fig = plot_top_features(deg_df, n=top_n,
                                    title=f"Top {top_n} DEGs by |Log2FC|")
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.warning("No significant DEGs found with current thresholds.")

    # -- DEG Table --
    with tabs[tab_idx]:
        tab_idx += 1
        st.subheader("Differentially Expressed Features Table")

        if is_dual:
            stage_choice = st.radio(
                "Select stage", ["Stage 1", "Stage 2"],
                horizontal=True, key="table_stage"
            )
            deg_df = results["degs_stage1"] if "Stage 1" in stage_choice else results["degs_stage2"]
        else:
            deg_df = results.get("degs", pd.DataFrame())

        if not deg_df.empty:
            display_df = deg_df.reset_index().rename(
                columns={"index": "feature_id", "feature_id": "feature_id"}
            )
            display_df["padj"] = display_df["padj"].apply(lambda x: f"{x:.2e}")
            display_df["pvalue"] = display_df["pvalue"].apply(lambda x: f"{x:.2e}")
            for col in ["baseMean", "log2FoldChange", "lfcSE", "stat"]:
                if col in display_df.columns:
                    display_df[col] = display_df[col].round(4)
            st.dataframe(display_df, use_container_width=True, height=450)

            csv = deg_df.reset_index().to_csv(index=False)
            st.download_button(
                "Download DEG Table (CSV)",
                data=csv,
                file_name=f"{st.session_state.project_id}_DEGs.csv",
                mime="text/csv",
            )
        else:
            st.warning("No significant DEGs found with current thresholds.")

    # -- Venn Diagram (dual stage only) --
    if is_dual:
        with tabs[tab_idx]:
            tab_idx += 1
            st.subheader("Gene Overlap: Tumor-associated vs. Clinically Susceptible")

            s1_genes = set(results["degs_stage1"].index)
            s2_genes = set(results["degs_stage2"].index)

            col1, col2, col3 = st.columns([1, 2, 1])
            with col2:
                try:
                    venn_png = plot_venn(
                        s1_genes, s2_genes,
                        label_a="Tumor-associated",
                        label_b="Clinically susceptible",
                        title="Shared DEGs",
                    )
                    st.image(venn_png, caption="Venn diagram of overlapping DEGs.")
                except ImportError as e:
                    st.warning(str(e))

            col1, col2, col3 = st.columns(3)
            col1.metric("Tumor-associated only", len(s1_genes - s2_genes))
            col2.metric("Shared", len(s1_genes & s2_genes))
            col3.metric("Clinically susceptible only", len(s2_genes - s1_genes))

        # -- Wald Scatter --
        with tabs[tab_idx]:
            tab_idx += 1
            st.subheader("Wald Test Value vs. Log2 Fold Change")

            res_s1 = results["results_stage1"]
            overlap_genes = (
                list(results["overlap"].index)
                if not results["overlap"].empty
                else []
            )
            fig = plot_wald_scatter(
                res_s1,
                highlight_genes=overlap_genes[:10],
                title="Stage 1: LFC vs. Wald Test (overlap genes highlighted)",
            )
            st.plotly_chart(fig, use_container_width=True)
            st.caption(
                "Blue points: genes in the final overlap table. "
                "Statistical significance is conveyed by Wald test value magnitude."
            )

        # -- Overlap Table --
        with tabs[tab_idx]:
            tab_idx += 1
            st.subheader(f"Biomarker Candidates (Overlap, |LFC| >= {overlap_lfc_threshold})")

            overlap_df = results.get("overlap", pd.DataFrame())
            if not overlap_df.empty:
                disp = overlap_df.reset_index()
                disp["padj"] = disp["padj"].apply(lambda x: f"{x:.2e}")
                for col in ["baseMean", "log2FoldChange", "lfcSE", "stat"]:
                    if col in disp.columns:
                        disp[col] = disp[col].round(4)
                st.dataframe(disp, use_container_width=True)

                csv_overlap = overlap_df.reset_index().to_csv(index=False)
                st.download_button(
                    "Download Overlap Table (CSV)",
                    data=csv_overlap,
                    file_name=f"{st.session_state.project_id}_overlap_biomarkers.csv",
                    mime="text/csv",
                )

                st.markdown(
                    """
                    These genes are differentially expressed in **both**:
                    - Tumor vs. normal tissue (Stage 1)
                    - Between the selected clinical groups (Stage 2)

                    They represent candidate biomarkers linking the clinical variable
                    to tumor development.
                    """
                )
            else:
                st.info(
                    "No genes meet the overlap criteria at the current thresholds. "
                    "Try lowering the LFC or padj thresholds."
                )

    # -- Sample Overview --
    with tabs[tab_idx]:
        tab_idx += 1
        st.subheader("Sample Overview")

        meta = pl.filtered_metadata if pl.filtered_metadata is not None else (
            pl.sample_info.join(pl.clinical_df, on="case_id", how="left")
            if pl.sample_info is not None and pl.clinical_df is not None
            else pl.sample_info
        )

        if meta is not None:
            for col, label in [
                ("sample_type", "Sample Type"),
                ("gender", "Gender"),
                ("race", "Race / Ethnicity"),
                ("tobacco_smoking_status", "Smoking Status"),
            ]:
                if col in meta.columns:
                    st.plotly_chart(
                        plot_sample_distribution(meta, col, f"{label} Distribution"),
                        use_container_width=True,
                    )
        else:
            st.info("No metadata available.")

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------

st.divider()
st.caption(
    "TCGA Cancer Epigenomics Pipeline | "
    "Department of Bioinformatics, i3L University, Jakarta, Indonesia | "
    "Data source: NCI Genomic Data Commons (https://portal.gdc.cancer.gov)"
)
