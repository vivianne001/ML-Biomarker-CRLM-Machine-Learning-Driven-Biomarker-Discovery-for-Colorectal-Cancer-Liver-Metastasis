# ML-Biomarker-CRLM-Machine-Learning-Driven-Biomarker-Discovery-for-Colorectal-Cancer-Liver-Metastasis
Programming language: R (>=4.4.1)
Project Overview
This repository implements a reproducible workflow for genomic, transcriptomic, and clinical data integration using machine learning algorithms to discover and validate biomarkers for colorectal cancer liver metastasis (CRLM).
Table of Contents
├── Building Models
│   ├── Baseline Test
│   │   ├── heatmap of baseline test in testing cohort.pdf
│   │   ├── heatmap of baseline test in validation cohort.pdf
│   │   ├── multi learners classification result.csv
│   │   └── visualization of baseline test results in validation and testing cohorts.R
│   ├── Feature Filtering
│   │   ├── Calculating the variance of input features.R
│   │   ├── Regularization filtering.R
│   │   ├── feature_filter_set_seed.Rdata
│   │   ├── filter-variance.Rdata
│   │   └── mlr3_classification_learners_list.Rdata
│   └── Model
│       ├── Candidate biomarkers selection and visualization.R
│       ├── Heatmap visuliation of  29 Candidate Biomarkers.R
│       ├── ML models for classification task.R
│       ├── Optimal K Value for KNN Model.R
│       ├── Visualization of Confusion Matrix.R
│       ├── acc_auc of testset.Rdata
│       ├── acc_auc of validationset.Rdata
│       ├── cohort split for model training and testing.Rdata
│       ├── confusion_matrix of validation datasets.Rdata
│       ├── feature_auc_values-permutation-trainset 29 + CRC BIOMARKER Features.Rdata
│       ├── feature_auc_values-testset validation-permutation-trainset 29.Rdata
│       ├── filter-variance.Rdata
│       └── removebatcheffects.trainset-all genes plus status.Rdata
├── Data Preprocessing
│   ├── Functional DEGs analysis
│   │   ├── GSE131418_Raw_Data.Rdata
│   │   ├── dataset preprocessing.R
│   │   └── group_metastasis_primary.Rdata
│   ├── Model datasets
│   │   ├── Cohort Split.R
│   │   ├── Dataset preprocessing.R
│   │   ├── GSE18105
│   │   │   ├── GSE18105_Raw_Data.Rdata
│   │   │   └── group_metastasis_primary.Rdata
│   │   ├── GSE21510
│   │   │   ├── GSE21510_Raw_Data.Rdata
│   │   │   └── group_metastasis_primary.Rdata
│   │   ├── GSE27854
│   │   │   ├── GSE27854_Raw_Data.Rdata
│   │   │   └── group_metastasis_primary.Rdata
│   │   ├── GSE41568
│   │   │   ├── GSE41568_Raw_Data.Rdata
│   │   │   └── group_metastasis_primary.Rdata
│   │   ├── GSE71222
│   │   │   ├── GSE71222_Raw_Data.Rdata
│   │   │   └── group_metastasis_primary.Rdata
│   │   ├── GSE81986
│   │   │   ├── GSE81986_Adjusted Data_Gene Symbol.Rdata
│   │   │   └── grouplist_metastasis_primary.Rdata
│   │   ├── Remove Batch Effect.R
│   │   ├── removebatcheffects.trainset.Rdata
│   │   └── trainset_validationset_testset.Rdata
│   └── Other datasets
│       ├── Dataset preprocessing.R
│       ├── GSE103479
│       │   ├── GPL23985-ID.xlsx
│       │   └── GSE103479_Raw_Data.Rdata
│       ├── GSE108277
│       │   ├── GSE108277_processed.Rdata
│       │   └── ID.xlsx
│       ├── GSE159216
│       │   ├── GSE159216 ID.xlsx
│       │   ├── GSE159216_Raw_Data.Rdata
│       │   ├── GSE159216_processed.Rdata
│       │   ├── exp matrix.csv
│       │   └── ph.csv
│       ├── GSE17536
│       │   ├── GSE17536_processed.Rdata.Rdata
│       │   └── GSE39582 ID.xlsx
│       ├── GSE204805
│       │   ├── GSE204805_merged_hs_mm.tsv.gz
│       │   ├── GSE204805_processed.Rdata
│       │   └── group_metastasis_primary.Rdata
│       ├── GSE28702
│       │   ├── GSE28702_processed.Rdata
│       │   └── GSE39582 ID.xlsx
│       ├── GSE50760
│       │   ├── GSE50760-EXP.xlsx
│       │   └── GSE50760_processed.Rdata
│       ├── GSE72970
│       │   ├── GSE39582 ID.xlsx
│       │   └── GSE72970_processed.Rdata
│       └── TCGA
│           ├── Dataset preprocessing.R
│           ├── TCGA.COADREAD.Survival time.xlsx
│           ├── TCGA.COADREAD.sampleMap_HiSeqV2.gz
│           └── adjusted raw data of CRC TCGA and grouplist.Rdata
└── Experiments
    ├── Biomarker Selection
    │   ├── GSE204805
    │   │   ├── Comparison of candidate biomarkers.R
    │   │   └── biomarker selection-Analysis of intergroup boxplots.R
    │   └── GSE50760
    │       ├── Comparison of candidate biomarkers.R
    │       └── biomarker selection-Analysis of intergroup boxplots.R
    ├── Clinical Analysis
    │   ├── Gene Expression Analysis Across Pathological Stages in TCGA.R
    │   ├── Kaplan-Meier Survival Analysis.R
    │   └── TGCA_COADREAD_Raw Counts.Rdata
    ├── Drug Therapeutic Prediction
    │   ├── Oncopredict
    │   │   ├── ACMSD group_high and low.Rdata
    │   │   ├── Drug Sensitivity Prediction Using OncoPredict.R
    │   │   ├── Oncopredict-GSE204805 BRAF drugs.Rdata
    │   │   ├── Oncopredict-TCGA_testPtype_NETS p38 CTRP2 drugs.Rdata
    │   │   ├── Training Datasets.Rdata
    │   │   ├── Visualization of Drug Sensitivity by OncoPredict.R
    │   │   ├── drug1 info of oncopredict training set.csv
    │   │   ├── oncopredict-drug and drug category.Rdata
    │   │   └── oncopredict_results.csv
    │   └── Other drug predictions
    │       ├── Drug Sensitivity analysis in GSE108277.R
    │       ├── Drug Sensitivity analysis in GSE28702.R
    │       └── Drug Sensitivity analysis in GSE72970.R
    ├── Functional DEGs
    │   ├── Differentially expression analysis and Pathway Enrichment.R
    │   ├── GSE131418-DEGS.csv
    │   ├── Identification and Visualization of Functional DEGs.R
    │   ├── function dataset.csv
    │   ├── functional dataset.Rdata
    │   ├── grouplist-volcano-metastasis-DATA.P.Value 0.05.Rdata
    │   └── metabolic genes.csv
    ├── Tissue Validation
    │   ├── 5 year post surgery recurrent or no recurrent.R
    │   ├── IHC Score quantitative analysis.R
    │   └── Masson Fibrosis quantitative analysis.R
    └── immune and fibrosis analysis
        ├── ACMSD group_high and low.Rdata
        ├── ACMSD high and low TIDE.csv
        ├── All_hg19gene_len.csv
        ├── Analysis of TIDE prediction scores.R
        ├── Barplot of Immune GSEA Pathways.R
        ├── Calculation of ESTIMATE Scores for Tumor Microenvironment Analysis.R
        ├── Correlation Analysis of Gene Expression with Tumor Purity.R
        ├── DEG of gsva_go_matrix.Rdata
        ├── GOPlot Analysis and Visualization for GSE204805.R
        ├── GOPlot circ of ACMSD high and low.rdata
        ├── GSEA analysis of both KEGG and GO database.R
        ├── GSVA Analysis and Visualization.R
        ├── GeneList.xlsx
        ├── Heatmap Analysis of Fibrosis and TGFB Pathway Markers.R
        ├── Immune_Inflammatory_GSEA_Statistics.xlsx
        ├── TGFB and Fibrosis markers.rdata
        ├── cor of ACMSD in CRLM samples.rdata
        ├── estimate.tpm.Rdata
        ├── exp of tpm.Rdata
        ├── gsva_go_matrix.csv
        ├── gsva_kegg_matrix.csv
        ├── msigdbr-KEGG-GO-pathways.Rdata
        ├── multibox plot analysis of inflammatory and immune markers.R
        ├── ssGSEA Analysis of the Tumor Immune Microenvironment.R
        └── ssGSEA_tpm_res.Rdata
  - 
