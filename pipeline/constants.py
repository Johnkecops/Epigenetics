"""
Constants: TCGA project codes, data type configurations, clinical field mappings.
"""

# All 33 TCGA cancer projects with full names
CANCER_TYPES = {
    "TCGA-ACC":  "Adrenocortical Carcinoma",
    "TCGA-BLCA": "Bladder Urothelial Carcinoma",
    "TCGA-BRCA": "Breast Invasive Carcinoma",
    "TCGA-CESC": "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma",
    "TCGA-CHOL": "Cholangiocarcinoma",
    "TCGA-COAD": "Colon Adenocarcinoma",
    "TCGA-DLBC": "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
    "TCGA-ESCA": "Esophageal Carcinoma",
    "TCGA-GBM":  "Glioblastoma Multiforme",
    "TCGA-HNSC": "Head and Neck Squamous Cell Carcinoma",
    "TCGA-KICH": "Kidney Chromophobe",
    "TCGA-KIRC": "Kidney Renal Clear Cell Carcinoma",
    "TCGA-KIRP": "Kidney Renal Papillary Cell Carcinoma",
    "TCGA-LAML": "Acute Myeloid Leukemia",
    "TCGA-LGG":  "Brain Lower Grade Glioma",
    "TCGA-LIHC": "Liver Hepatocellular Carcinoma",
    "TCGA-LUAD": "Lung Adenocarcinoma",
    "TCGA-LUSC": "Lung Squamous Cell Carcinoma",
    "TCGA-MESO": "Mesothelioma",
    "TCGA-OV":   "Ovarian Serous Cystadenocarcinoma",
    "TCGA-PAAD": "Pancreatic Adenocarcinoma",
    "TCGA-PCPG": "Pheochromocytoma and Paraganglioma",
    "TCGA-PRAD": "Prostate Adenocarcinoma",
    "TCGA-READ": "Rectum Adenocarcinoma",
    "TCGA-SARC": "Sarcoma",
    "TCGA-SKCM": "Skin Cutaneous Melanoma",
    "TCGA-STAD": "Stomach Adenocarcinoma",
    "TCGA-TGCT": "Testicular Germ Cell Tumors",
    "TCGA-THCA": "Thyroid Carcinoma",
    "TCGA-THYM": "Thymoma",
    "TCGA-UCEC": "Uterine Corpus Endometrial Carcinoma",
    "TCGA-UCS":  "Uterine Carcinosarcoma",
    "TCGA-UVM":  "Uveal Melanoma",
}

# GDC API data type configurations
DATA_TYPE_CONFIG = {
    "Gene Expression": {
        "data_category": "Transcriptome Profiling",
        "data_type": "Gene Expression Quantification",
        "workflow_type": "STAR - Counts",
        "count_column": "unstranded",
        "id_column": "gene_id",
        "name_column": "gene_name",
        "skip_rows": 1,       # Skip the header summary rows in STAR output
        "summary_prefixes": ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"],
    },
    "miRNA Expression": {
        "data_category": "Transcriptome Profiling",
        "data_type": "miRNA Expression Quantification",
        "workflow_type": None,
        "count_column": "read_count",
        "id_column": "miRNA_ID",
        "name_column": "miRNA_ID",
        "skip_rows": 0,
        "summary_prefixes": [],
    },
}

# GDC API base URL
GDC_BASE_URL = "https://api.gdc.cancer.gov"

# TCGA barcode sample type codes
# Codes 01-09 = tumor; 10-19 = normal solid tissue; 20+ = control/cell line
TUMOR_CODES = {str(i).zfill(2) for i in range(1, 10)}
NORMAL_CODES = {str(i).zfill(2) for i in range(10, 20)}

# GDC clinical field paths for metadata
CLINICAL_FIELDS = [
    "case_id",
    "submitter_id",
    "demographic.gender",
    "demographic.race",
    "demographic.ethnicity",
    "demographic.vital_status",
    "demographic.age_at_index",
    "exposures.tobacco_smoking_status",
    "exposures.tobacco_smoking_history_indicator",
    "exposures.pack_years_smoked",
    "diagnoses.primary_diagnosis",
    "diagnoses.ajcc_pathologic_stage",
    "diagnoses.tumor_stage",
    "samples.sample_type",
    "samples.sample_type_id",
    "samples.submitter_id",
]

# Human-readable smoking status labels (TCGA tobacco_smoking_status codes)
SMOKING_STATUS_LABELS = {
    "1": "Lifelong Non-Smoker (<100 cigarettes)",
    "2": "Current Smoker",
    "3": "Current Reformed Smoker (>15 years)",
    "4": "Current Reformed Smoker (<=15 years)",
    "5": "Current Reformed Smoker (Duration Unknown)",
    "6": "Smoker at Diagnosis",
    "7": "Smoking History Not Documented",
    "Unknown": "Unknown",
}

# DESeq2 default parameters (mirrors paper methodology)
DEFAULT_PARAMS = {
    "lfc_threshold": 1.0,        # log2 fold change threshold for DEG calling
    "padj_threshold": 0.01,      # adjusted p-value threshold
    "min_count": 10,             # minimum count sum across all samples
    "min_samples_fraction": 0.5, # fraction of samples that must express gene
    "overlap_lfc_threshold": 2.0,# LFC threshold for final overlap table (paper Table 2)
    "top_n_genes": 5,            # top N genes to highlight (paper reports top 5)
}

# API pagination
GDC_PAGE_SIZE = 100
GDC_MAX_RETRIES = 3
GDC_RETRY_DELAY = 5  # seconds
