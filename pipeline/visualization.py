"""
visualization.py
----------------
Publication-quality interactive plots for the TCGA epigenomics pipeline.

Implements all visualizations from the reference papers:
- Figure 1/3 from Parikesit et al. (2023): MA plot (mean normalized counts vs LFC)
- Figure 2/4 from Parikesit et al. (2023): Volcano plot (LFC vs -log10 padj)
- Figure 5 from Parikesit et al. (2023): Venn diagram
- Figure 6 from Parikesit et al. (2023): Scatter plot (LFC vs Wald stat)
- Figure 1 from Fugaha et al. (2024): Volcano plot for miRNA

Uses Plotly for interactive Streamlit rendering; Matplotlib for Venn diagrams.
"""

import io
import logging
from typing import List, Optional, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

logger = logging.getLogger(__name__)

# Color scheme consistent across papers (blue = significant, grey = unchanged)
COLOR_ABERRANT = "#2166ac"
COLOR_UNCHANGED = "#b2b2b2"
COLOR_UP = "#d73027"
COLOR_DOWN = "#4575b4"


# ---------------------------------------------------------------------------
# MA plot
# ---------------------------------------------------------------------------

def plot_ma(
    results: pd.DataFrame,
    lfc_threshold: float = 1.0,
    padj_threshold: float = 0.01,
    highlight_genes: Optional[List[str]] = None,
    title: str = "MA Plot",
) -> go.Figure:
    """
    MA plot: mean normalized counts (x) vs. log2 fold change (y).

    Significant genes colored blue; threshold lines shown.
    Mirrors Figure 1 and Figure 3 from Parikesit et al. (2023).
    """
    df = results.copy().reset_index()
    df = df.dropna(subset=["baseMean", "log2FoldChange", "padj"])
    df["log_mean"] = np.log10(df["baseMean"].clip(lower=0.01))

    sig_mask = (df["padj"] < padj_threshold) & (df["log2FoldChange"].abs() >= lfc_threshold)
    df["color_group"] = np.where(sig_mask, "Aberrant", "Unchanged")

    id_col = df.columns[0]  # feature_id or miRNA_ID

    fig = px.scatter(
        df,
        x="log_mean",
        y="log2FoldChange",
        color="color_group",
        color_discrete_map={"Aberrant": COLOR_ABERRANT, "Unchanged": COLOR_UNCHANGED},
        hover_data={id_col: True, "baseMean": ":.2f", "padj": ":.2e",
                    "log_mean": False, "color_group": False},
        opacity=0.6,
        title=title,
        labels={
            "log_mean": "log10(Mean Normalized Counts)",
            "log2FoldChange": "Log Fold Change",
            "color_group": "Expression",
        },
    )

    fig.add_hline(y=0, line_dash="dash", line_color="black", line_width=1)

    # Highlight specific genes
    if highlight_genes:
        hl = df[df[id_col].isin(highlight_genes)]
        fig.add_scatter(
            x=hl["log_mean"], y=hl["log2FoldChange"],
            mode="markers+text",
            marker=dict(color=COLOR_UP, size=10, symbol="diamond"),
            text=hl[id_col],
            textposition="top right",
            name="Highlighted",
        )

    fig.update_layout(
        template="plotly_white",
        legend_title_text="Expression",
        height=450,
        font=dict(family="Arial", size=12),
    )
    return fig


# ---------------------------------------------------------------------------
# Volcano plot
# ---------------------------------------------------------------------------

def plot_volcano(
    results: pd.DataFrame,
    lfc_threshold: float = 1.0,
    padj_threshold: float = 0.01,
    highlight_genes: Optional[List[str]] = None,
    title: str = "Volcano Plot",
) -> go.Figure:
    """
    Volcano plot: log2 fold change (x) vs. -log10(adjusted p-value) (y).

    Vertical lines at ±lfc_threshold; horizontal line at -log10(padj_threshold).
    Mirrors Figures 2 and 4 from Parikesit et al. (2023) and Figure 1 from
    Fugaha et al. (2024).
    """
    df = results.copy().reset_index()
    df = df.dropna(subset=["log2FoldChange", "padj"])
    df["padj_clipped"] = df["padj"].clip(lower=1e-300)
    df["neg_log10_padj"] = -np.log10(df["padj_clipped"])

    sig_mask = (df["padj"] < padj_threshold) & (df["log2FoldChange"].abs() >= lfc_threshold)
    df["color_group"] = "Unchanged"
    df.loc[sig_mask & (df["log2FoldChange"] > 0), "color_group"] = "Overexpressed"
    df.loc[sig_mask & (df["log2FoldChange"] < 0), "color_group"] = "Underexpressed"

    id_col = df.columns[0]

    color_map = {
        "Overexpressed": COLOR_UP,
        "Underexpressed": COLOR_DOWN,
        "Unchanged": COLOR_UNCHANGED,
    }

    fig = px.scatter(
        df,
        x="log2FoldChange",
        y="neg_log10_padj",
        color="color_group",
        color_discrete_map=color_map,
        hover_data={id_col: True, "log2FoldChange": ":.3f",
                    "padj": ":.2e", "neg_log10_padj": False, "color_group": False},
        opacity=0.6,
        title=title,
        labels={
            "log2FoldChange": "Log2 Fold Change",
            "neg_log10_padj": "-log10(Adjusted P-value)",
            "color_group": "Expression",
        },
        category_orders={"color_group": ["Overexpressed", "Underexpressed", "Unchanged"]},
    )

    # Threshold lines
    fig.add_vline(x=lfc_threshold, line_dash="dash", line_color=COLOR_UP, line_width=1.5)
    fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color=COLOR_DOWN, line_width=1.5)
    fig.add_hline(
        y=-np.log10(padj_threshold),
        line_dash="dash", line_color="red", line_width=1.5
    )

    # Highlight specific genes
    if highlight_genes:
        hl = df[df[id_col].isin(highlight_genes)]
        fig.add_scatter(
            x=hl["log2FoldChange"], y=hl["neg_log10_padj"],
            mode="markers+text",
            marker=dict(color="orange", size=12, symbol="diamond"),
            text=hl[id_col],
            textposition="top right",
            name="Highlighted",
        )

    fig.update_layout(
        template="plotly_white",
        height=450,
        font=dict(family="Arial", size=12),
    )
    return fig


# ---------------------------------------------------------------------------
# Wald stat scatter plot
# ---------------------------------------------------------------------------

def plot_wald_scatter(
    results: pd.DataFrame,
    highlight_genes: Optional[List[str]] = None,
    title: str = "Log2 Fold Change vs. Wald Test Value",
) -> go.Figure:
    """
    Scatter plot of log2 fold change (x) vs. Wald test statistic (y).

    Highlights top implicated genes.
    Mirrors Figure 6 from Parikesit et al. (2023).
    """
    df = results.copy().reset_index()
    df = df.dropna(subset=["log2FoldChange", "stat"])

    id_col = df.columns[0]
    df["highlighted"] = df[id_col].isin(highlight_genes or [])

    fig = go.Figure()

    # Background points
    bg = df[~df["highlighted"]]
    fig.add_trace(go.Scatter(
        x=bg["log2FoldChange"],
        y=bg["stat"],
        mode="markers",
        marker=dict(color=COLOR_UNCHANGED, size=4, opacity=0.5),
        name="Not highlighted",
        text=bg[id_col],
        hovertemplate="%{text}<br>LFC: %{x:.3f}<br>Wald: %{y:.3f}<extra></extra>",
    ))

    # Highlighted points
    if highlight_genes:
        hl = df[df["highlighted"]]
        fig.add_trace(go.Scatter(
            x=hl["log2FoldChange"],
            y=hl["stat"],
            mode="markers+text",
            marker=dict(color=COLOR_ABERRANT, size=10, symbol="circle"),
            text=hl[id_col],
            textposition="top right",
            name="Implicated genes",
            hovertemplate="%{text}<br>LFC: %{x:.3f}<br>Wald: %{y:.3f}<extra></extra>",
        ))

    fig.add_hline(y=0, line_dash="dash", line_color="black", line_width=1)
    fig.add_vline(x=0, line_dash="dash", line_color="black", line_width=1)

    fig.update_layout(
        title=title,
        xaxis_title="Log2 Fold Change",
        yaxis_title="Wald Test Value",
        template="plotly_white",
        height=450,
        font=dict(family="Arial", size=12),
    )
    return fig


# ---------------------------------------------------------------------------
# Venn diagram
# ---------------------------------------------------------------------------

def plot_venn(
    set_a: Set[str],
    set_b: Set[str],
    label_a: str = "Group A",
    label_b: str = "Group B",
    title: str = "Gene Overlap",
) -> bytes:
    """
    Two-set Venn diagram using matplotlib_venn.

    Returns PNG bytes for display in Streamlit.
    Mirrors Figure 5 from Parikesit et al. (2023).
    """
    try:
        from matplotlib_venn import venn2
    except ImportError:
        raise ImportError(
            "matplotlib_venn is required for Venn diagrams. "
            "Install with: pip install matplotlib-venn"
        )

    fig, ax = plt.subplots(figsize=(6, 5))
    v = venn2(
        [set_a, set_b],
        set_labels=(label_a, label_b),
        ax=ax,
        set_colors=("#4575b4", "#d73027"),
        alpha=0.5,
    )

    if v.get_label_by_id("10"):
        v.get_label_by_id("10").set_fontsize(14)
    if v.get_label_by_id("11"):
        v.get_label_by_id("11").set_fontsize(14)
        v.get_label_by_id("11").set_fontweight("bold")
    if v.get_label_by_id("01"):
        v.get_label_by_id("01").set_fontsize(14)

    ax.set_title(title, fontsize=14, pad=12)
    plt.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf.read()


# ---------------------------------------------------------------------------
# Sample distribution bar chart
# ---------------------------------------------------------------------------

def plot_sample_distribution(
    metadata: pd.DataFrame,
    group_col: str = "sample_type",
    title: str = "Sample Distribution",
) -> go.Figure:
    """Bar chart of sample counts grouped by a metadata column."""
    if group_col not in metadata.columns:
        fig = go.Figure()
        fig.update_layout(title=f"Column '{group_col}' not available.")
        return fig

    counts = metadata[group_col].value_counts().reset_index()
    counts.columns = [group_col, "count"]

    fig = px.bar(
        counts,
        x=group_col,
        y="count",
        color=group_col,
        title=title,
        labels={group_col: group_col.replace("_", " ").title(), "count": "Sample Count"},
        text="count",
    )
    fig.update_traces(textposition="outside")
    fig.update_layout(
        template="plotly_white",
        showlegend=False,
        height=350,
        font=dict(family="Arial", size=12),
    )
    return fig


# ---------------------------------------------------------------------------
# Top features bar chart
# ---------------------------------------------------------------------------

def plot_top_features(
    degs: pd.DataFrame,
    n: int = 20,
    value_col: str = "log2FoldChange",
    title: str = "Top Differentially Expressed Features",
) -> go.Figure:
    """
    Horizontal bar chart of top N features sorted by LFC magnitude.
    """
    df = degs.copy().reset_index()
    id_col = df.columns[0]

    df["abs_lfc"] = df[value_col].abs()
    df = df.sort_values("abs_lfc", ascending=False).head(n)
    df = df.sort_values(value_col)

    df["color"] = np.where(df[value_col] > 0, COLOR_UP, COLOR_DOWN)

    fig = go.Figure(go.Bar(
        x=df[value_col],
        y=df[id_col],
        orientation="h",
        marker_color=df["color"],
        text=df[value_col].round(2),
        textposition="outside",
        hovertemplate="%{y}<br>LFC: %{x:.3f}<extra></extra>",
    ))

    fig.add_vline(x=0, line_color="black", line_width=1)
    fig.update_layout(
        title=title,
        xaxis_title="Log2 Fold Change",
        yaxis_title="",
        template="plotly_white",
        height=max(300, n * 22),
        font=dict(family="Arial", size=11),
    )
    return fig


# ---------------------------------------------------------------------------
# Heatmap of top features
# ---------------------------------------------------------------------------

def plot_expression_heatmap(
    count_matrix: pd.DataFrame,
    gene_list: List[str],
    metadata: pd.DataFrame,
    group_col: str = "condition",
    title: str = "Expression Heatmap (Top DEGs)",
    max_samples: int = 100,
) -> go.Figure:
    """
    Z-score normalized heatmap of top DEGs across samples.
    Samples sorted by group label.
    """
    genes_present = [g for g in gene_list if g in count_matrix.columns]
    if not genes_present:
        fig = go.Figure()
        fig.update_layout(title="No matching genes in count matrix.")
        return fig

    meta = metadata.copy()
    if group_col in meta.columns:
        meta = meta.sort_values(group_col)

    samples = meta.index.intersection(count_matrix.index)
    if len(samples) > max_samples:
        samples = samples[:max_samples]

    mat = count_matrix.loc[samples, genes_present].T.astype(float)

    # Z-score per gene
    mat = mat.subtract(mat.mean(axis=1), axis=0)
    std = mat.std(axis=1).replace(0, 1)
    mat = mat.divide(std, axis=0)

    group_labels = (
        meta.loc[samples, group_col].tolist() if group_col in meta.columns else None
    )

    fig = go.Figure(go.Heatmap(
        z=mat.values,
        x=samples.tolist(),
        y=genes_present,
        colorscale="RdBu_r",
        zmid=0,
        colorbar=dict(title="Z-score"),
        hovertemplate="Gene: %{y}<br>Sample: %{x}<br>Z-score: %{z:.2f}<extra></extra>",
    ))

    fig.update_layout(
        title=title,
        xaxis_title="Samples",
        yaxis_title="Features",
        template="plotly_white",
        height=max(400, len(genes_present) * 25),
        font=dict(family="Arial", size=11),
    )
    return fig
