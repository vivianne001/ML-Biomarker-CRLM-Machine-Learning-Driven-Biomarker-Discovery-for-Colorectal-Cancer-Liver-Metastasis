# ML-Biomarker-CRLM-Machine-Learning-Driven-Biomarker-Discovery-for-Colorectal-Cancer-Liver-Metastasis
Programming language: R (>=4.4.1)
Project Overview
This repository implements a reproducible workflow for genomic, transcriptomic, and clinical data integration using machine learning algorithms to discover and validate biomarkers for colorectal cancer liver metastasis (CRLM).
________________________________________________________
Table of Contents
├── Building Models
│   ├── Baseline Test
│   ├── Feature Filtering
│   └── Model
├── Data Preprocessing
│   ├── Functional DEGs analysis
│   ├── Model datasets
│   └── Other datasets
└── Experiments
    ├── Biomarker Selection
    ├── Clinical Analysis
    ├── Drug Therapeutic Prediction
    ├── Functional DEGs
    ├── Tissue Validation
    └── immune and fibrosis analysis
___________________________________________________
Key Features & Highlights
- Large-cohort integration: Harmonize >1,600 CRC samples (TCGA, GEO, plus original tissues samples).
- Flexible machine learning modules: Fast feature screening and classification using 20 algorithms (mlr3, caret, etc.).
- Comprehensive biological/clinical annotation: Downstream enrichment, immune microenvironment, survival, drug response, and fibrosis analysis.
- Reproducibility: Full code, statistics, and visualization environments are documented.
- Extensible: Modular design allows adaptation to other cancer biomarker and risk modeling projects.

## Installation & Dependencies
Required:

- R >= 4.4.1
- Core R packages: [mlr3, mlr3extralearners, limma, clusterProfiler, GSVA, survminer, pROC, estimate, ComplexHeatmap, oncoPredict, MetaboAnalystR, etc.]

Install dependencies:
# CRAN
install.packages(c("mlr3", "mlr3extralearners", "pROC", "survminer", "ggplot2", "data.table", ...))
# Bioconductor
if(!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("limma", "GSVA", "ComplexHeatmap", "clusterProfiler", "estimate"))

Full dependency list and version info is provided via sessionInfo.txt in the repository.

# Data Sources
- Public transcriptomic datasets (TCGA, GEO: accession numbers provided in Supplementary Table S1,from manuscript)
- In-house tissue/serum validation data (details and access: see manuscript Methods/Data Availability)
- Clinical and metabolomics metadata (see manuscript for details)

Note: Data use is for academic and non-commercial research only, in accordance with original database licenses.

