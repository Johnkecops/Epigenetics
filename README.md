# TCGA Cancer Epigenomics Pipeline

A Python bioinformatics pipeline for differential gene and miRNA expression analysis using data from the NCI Cancer Genome Atlas (TCGA), with an interactive Streamlit front-end.

## Overview

This tool implements the two-stage differential expression methodology described in:

- Parikesit AA et al. *Comparative Gene Analysis of Squamous Cell Lung Carcinoma Between Smoking and Non-smoking Individuals*. Berita Biologi 22(3), 2023. DOI: 10.14203/beritabiology.v20i1.3991
- Fugaha DR, Putta D, Lunoto D, Parikesit AA. *MicroRNA Expression Profiling in Colorectal Carcinoma: Identification and Validation of Novel Diagnostic, Prognostic and Predictive Biomarkers*. ICIC 2024.

The pipeline retrieves count data and clinical metadata directly from the GDC API, applies user-selected clinical filters, runs DESeq2 differential expression analysis via `pydeseq2`, and generates interactive visualizations.

## Features

- **Universal cancer support**: All 33 TCGA cancer project codes (LUSC, BRCA, COAD, LIHC, PRAD, and 28 others)
- **Two molecular layers**: Gene expression (mRNA, STAR counts) and miRNA expression
- **Clinical stratification**: Filter patients by gender, race/ethnicity, and smoking status before analysis
- **Two analysis modes**:
  - *Single-stage*: Tumor vs. normal (Fugaha et al. 2024 style)
  - *Dual-stage*: Tumor/normal + clinical group comparison with gene overlap (Parikesit et al. 2023 style)
- **Visualizations**: Volcano plot, MA plot, Wald scatter plot, Venn diagram, expression heatmap
- **Streamlit front-end**: Interactive, step-by-step web interface

## Installation

### 1. Create a virtual environment

```bash
python -m venv venv
source venv/bin/activate        # macOS / Linux
venv\Scripts\activate.bat       # Windows
```

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

> **R is not required.** This pipeline uses `pydeseq2`, a pure-Python port of the R DESeq2 package.

## Usage

### Streamlit App (recommended)

```bash
streamlit run app.py
```

Open your browser at `http://localhost:8501`.

**Workflow:**
1. Select cancer type (e.g., TCGA-LUSC) and data layer (Gene Expression / miRNA)
2. Click **Load Data from TCGA** - data is downloaded from GDC and cached locally
3. Apply clinical filters: gender, race/ethnicity, smoking status
4. Click **Apply Filters**
5. Choose analysis mode and set DESeq2 thresholds
6. Click **Run Analysis**
7. Explore volcano plot, MA plot, DEG tables, Venn diagram, and overlap biomarkers
8. Download results as CSV

### Python API

```python
from pipeline.pipeline import EpigenomicsPipeline

# Initialize with local cache directory
pl = EpigenomicsPipeline(cache_dir="./tcga_cache")

# Step 1: Download data
pl.load_data("TCGA-LUSC", "Gene Expression")

# Step 2: Apply clinical filters
pl.set_clinical_filters({
    "gender": ["male"],
    "tobacco_smoking_status": ["2", "3", "4"],  # smokers
})

# Step 3: Run dual-stage analysis (Parikesit et al. 2023 methodology)
results = pl.run_analysis(
    analysis_mode="dual_stage",
    stage2_col="tobacco_smoking_status",
    stage2_group_a=["2", "3", "4"],   # smokers
    stage2_group_b=["1"],              # non-smokers
    stage2_label_a="Smoker",
    stage2_label_b="Non-smoker",
    lfc_threshold=1.0,
    padj_threshold=0.01,
    overlap_lfc_threshold=2.0,
)

# Results
print(results["summary"])
print(results["overlap"])          # Biomarker candidates
print(results["degs_stage1"])      # Stage 1 DEGs
```

## Pipeline Architecture

```
app.py                          Streamlit front-end
pipeline/
    constants.py                Cancer types, data configs, field mappings
    tcga_downloader.py          GDC REST API wrapper (file query + download + parse)
    preprocessing.py            Clinical filters, tumor/normal labeling, expression filtering
    deseq_analysis.py           pydeseq2 wrapper, dual/single-stage orchestration
    visualization.py            Plotly interactive plots + matplotlib Venn diagrams
    pipeline.py                 EpigenomicsPipeline class (orchestrates all modules)
```

### Analysis workflow (dual-stage)

```
TCGA GDC API
    -> query_files()            [query file catalog by project + data type]
    -> get_clinical_metadata()  [gender, race, smoking status, tumor stage]
    -> download_and_parse()     [batch download + count matrix assembly]
         |
    apply_clinical_filters()    [restrict to subpopulation of interest]
         |
    Stage 1: DESeq2             [tumor vs. normal within clinical group]
         |                       design = ~ condition
    Stage 2: DESeq2             [clinical variable comparison]
         |                       design = ~ condition + tumor_status
    Overlap                     [intersect Stage 1 and Stage 2 DEG lists]
         |                       filter by |LFC| >= overlap_lfc_threshold
    Visualize                   [Volcano, MA, Wald scatter, Venn, heatmap]
```

## Data Source

All data is retrieved from the [NCI Genomic Data Commons (GDC)](https://portal.gdc.cancer.gov/) via the public REST API. No authentication is required for open-access TCGA data.

Downloaded files are cached locally under `./tcga_cache/` to avoid repeated downloads.

## DESeq2 Filtering Criteria

Following Parikesit et al. (2023) and Yao et al. (2019):

| Criterion | Default value |
|-----------|--------------|
| Minimum count per sample | 10 |
| Fraction of samples expressing | >= 50% |
| log2 fold change (DEG call) | >= 1.0 |
| Adjusted p-value (FDR) | < 0.01 |
| Overlap LFC threshold | >= 2.0 |

All thresholds are adjustable in the Streamlit sidebar.

## TCGA Smoking Status Codes

| Code | Meaning |
|------|---------|
| 1 | Lifelong Non-Smoker (<100 cigarettes) |
| 2 | Current Smoker |
| 3 | Current Reformed Smoker (>15 years) |
| 4 | Current Reformed Smoker (<=15 years) |
| 5 | Current Reformed Smoker (Duration Unknown) |

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@article{parikesit2023lusc,
  author  = {Shemuel, Josia and Angelique, Priscilla and Josephine, Evelina
             and Fugaha, Daniel Ryan and Gabriela, Vania
             and Alhusain, Shaheer and Parikesit, Arli Aditya},
  title   = {Comparative Gene Analysis of Squamous Cell Lung Carcinoma
             Between Smoking and Non-smoking Individuals},
  journal = {Berita Biologi},
  volume  = {22},
  number  = {3},
  pages   = {291--302},
  year    = {2023},
  doi     = {10.14203/beritabiology.v20i1.3991}
}

@inproceedings{fugaha2024mirna,
  author    = {Fugaha, Daniel Ryan and Putta, Dhannyo and Lunoto, Dennis
               and Parikesit, Arli Aditya},
  title     = {MicroRNA Expression Profiling in Colorectal Carcinoma:
               Identification and Validation of Novel Diagnostic,
               Prognostic and Predictive Biomarkers},
  booktitle = {2024 Ninth International Conference on Informatics and Computing (ICIC)},
  year      = {2024}
}
```

## License

MIT License. See [LICENSE.txt](LICENSE.txt).

## Author

**Dr.rer.nat. Arli Aditya Parikesit**
Department of Bioinformatics, Indonesia International Institute for Life Sciences (i3L)
Jakarta, Indonesia
ORCID: [0000-0001-8716-3926](https://orcid.org/0000-0001-8716-3926)
Email: arli.parikesit@i3l.ac.id
