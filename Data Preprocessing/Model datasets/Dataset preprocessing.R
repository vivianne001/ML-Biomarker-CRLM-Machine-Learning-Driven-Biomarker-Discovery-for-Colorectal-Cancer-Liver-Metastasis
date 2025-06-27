# ============================================================================
# GEO Datasets Preprocessing Script (GSE131418, GSE18105, GSE21510, GSE27854, GSE41568, GSE71222, GSE81986)
# ============================================================================

rm(list = ls())

# ===============================
# Required Libraries
# ===============================
library(BiocManager)
library(GEOquery)
library(limma)
library(readxl)
options(stringsAsFactors = FALSE)

# ===============================
# Helper Functions
# ===============================

normalize_expression <- function(exp_matrix) {
  boxplot(exp_matrix) # Check signal distribution
  print(range(exp_matrix))
  exp_matrix = normalizeBetweenArrays(exp_matrix)
  print(range(exp_matrix))
  # Log transformation if needed
  qx <- as.numeric(quantile(exp_matrix, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
  LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) {
    exp_matrix[which(exp_matrix <= 0)] <- NaN
    exp_matrix <- log2(exp_matrix+1)
  }
  print(range(exp_matrix))
  return(as.data.frame(exp_matrix))
}

map_gene_symbols <- function(exp, id_table, gene_col = "Gene Symbol") {
  id_table <- id_table[!is.na(id_table[[gene_col]]), ]
  exp$ID <- rownames(exp)
  exp <- merge(exp, id_table, by = "ID")
  exp <- as.data.frame(exp)
  # Remove ID and extra columns, set gene symbols as row names
  exp <- exp[, !(colnames(exp) %in% c("ID"))]
  rownames(exp) <- exp[[gene_col]]
  # Remove duplicated gene symbols
  exp <- exp[!duplicated(rownames(exp)), ]
  # Remove gene symbol column if present
  if (gene_col %in% colnames(exp)) {
    exp <- exp[, !(colnames(exp) %in% gene_col)]
  }
  return(exp)
}

set_groups <- function(ph, col, metastatic_value, group_names = c("metastasis", "primary")) {
  grouplist <- ifelse(ph[[col]] == metastatic_value, group_names[1], group_names[2])
  grouplist <- factor(grouplist, levels = group_names)
  return(grouplist)
}

filter_samples <- function(ph, exp, sample_id_col) {
  common_samples <- intersect(colnames(exp), ph[[sample_id_col]])
  ph <- ph[ph[[sample_id_col]] %in% common_samples, ]
  exp <- exp[, common_samples]
  if (!identical(ph[[sample_id_col]], colnames(exp))) {
    idx <- match(ph[[sample_id_col]], colnames(exp))
    exp <- exp[, idx]
  }
  return(list(ph = ph, exp = exp))
}

# ===============================
# Dataset-specific Preprocessing
# ===============================

# ---- GSE131418 ----
# Set working directory and read files as needed
setwd("path_to/GSE131418")
eSet <- getGEO(filename = "GSE131418_series_matrix.txt.gz", getGPL = FALSE)
data1 <- read.table("GSE131418_Consortium_prim_met_GE_probe_level.txt.gz", header = TRUE, sep = "\t")
data2 <- read.table("GSE131418_MCC_prim_met_GE_probe_level.txt.gz", header = TRUE, sep = "\t")
exp <- cbind(data1, data2)
exp <- normalize_expression(exp)

ph <- as.data.frame(pData(eSet))
ID <- as.data.frame(read_excel("GSE 131418 symbol.xlsx", sheet = 1, col_names = TRUE))
exp <- map_gene_symbols(exp, ID, gene_col = "GeneSymbol")
ph <- ph[ph$`treatment status classification (see description):ch1` %in% "PRE", ]
tmp <- filter_samples(ph, exp, "title")
ph <- tmp$ph; exp <- tmp$exp
grouplist <- ifelse(ph$characteristics_ch1.6 == "tumor type: METASTASIS", "metastasis", "primary")
grouplist <- factor(grouplist, levels = c("metastasis", "primary"))

save(exp, ID, ph, file = "GSE131418_Raw_Data.Rdata")
save(grouplist, file = "group_metastasis_primary.Rdata")

# ---- GSE18105 ----
setwd("path_to/GSE18105")
eSet <- getGEO(filename = "GSE18105_series_matrix.txt.gz", getGPL = FALSE)
exp <- exprs(eSet)
exp <- normalize_expression(exp)

ph <- as.data.frame(pData(eSet))
ID <- as.data.frame(read_excel("GPL570.xlsx", sheet = 1, col_names = TRUE))
exp <- map_gene_symbols(exp, ID, gene_col = "Gene Symbol")
ph <- ph[ph$`tissue:ch1` %in% c("cancer, homogenized ", "cancer, LCM"), ]
tmp <- filter_samples(ph, exp, "rownames(ph)")
ph <- tmp$ph; exp <- tmp$exp
# Remove anything after ' /// ' in row names and make unique
rownames(exp) <- make.unique(sub(" /// .*", "", rownames(exp)))
grouplist <- ifelse(ph$`metastasis:ch1` == "none", "primary", "metastasis")
grouplist <- factor(grouplist, levels = c("metastasis", "primary"))

save(exp, ID, ph, file = "GSE18105_Raw_Data.Rdata")
save(grouplist, file = "group_metastasis_primary.Rdata")

# ---- GSE21510 ----
setwd("path_to/GSE21510")
eSet <- getGEO(filename = "GSE21510_series_matrix.txt.gz", getGPL = FALSE)
exp <- exprs(eSet)
exp <- normalize_expression(exp)

ph <- as.data.frame(pData(eSet))
ID <- as.data.frame(read_excel("GPL570.xlsx", sheet = 1, col_names = TRUE))
exp <- map_gene_symbols(exp, ID, "Gene Symbol")
ph <- ph[ph$`tissue:ch1` %in% c("cancer, homogenized ", "cancer, LCM"), ]
tmp <- filter_samples(ph, exp, "rownames(ph)")
ph <- tmp$ph; exp <- tmp$exp
rownames(exp) <- make.unique(sub(" /// .*", "", rownames(exp)))
grouplist <- ifelse(ph$`metastatic tumor site:ch1` == "NA", "primary", "metastasis")
grouplist <- factor(grouplist, levels = c("metastasis", "primary"))

save(exp, ID, ph, file = "GSE21510_Raw_Data.Rdata")
save(grouplist, file = "group_metastasis_primary.Rdata")

# ---- GSE27854 ----
setwd("path_to/GSE27854")
eSet <- getGEO(filename = "GSE27854_series_matrix.txt.gz", getGPL = FALSE)
exp <- exprs(eSet)
exp <- normalize_expression(exp)

ph <- as.data.frame(pData(eSet))
ID <- as.data.frame(read_excel("GPL570.xlsx", sheet = 1, col_names = TRUE))
exp <- map_gene_symbols(exp, ID, "Gene Symbol")
ph <- ph[ph$`tissue:ch1` %in% c("cancer, homogenized ", "cancer, LCM"), ]
tmp <- filter_samples(ph, exp, "rownames(ph)")
ph <- tmp$ph; exp <- tmp$exp
rownames(exp) <- make.unique(sub(" /// .*", "", rownames(exp)))
grouplist <- ifelse(ph$`metastatic tumor site:ch1` == "NA", "primary", "metastasis")
grouplist <- factor(grouplist, levels = c("metastasis", "primary"))

save(exp, ID, ph, file = "GSE27854_Raw_Data.Rdata")
save(grouplist, file = "group_metastasis_primary.Rdata")

# ---- GSE41568 ----
setwd("path_to/GSE41568")
eSet <- getGEO(filename = "GSE41568_series_matrix.txt.gz", getGPL = FALSE)
exp <- exprs(eSet)
exp <- normalize_expression(exp)

ph <- as.data.frame(pData(eSet))
ID <- as.data.frame(read_excel("GPL570.xlsx", sheet = 1, col_names = TRUE))
exp <- map_gene_symbols(exp, ID, "Gene Symbol")
tmp <- filter_samples(ph, exp, "rownames(ph)")
ph <- tmp$ph; exp <- tmp$exp
rownames(exp) <- make.unique(sub(" /// .*", "", rownames(exp)))
grouplist <- ifelse(ph$`metastatic tumor site:ch1` == "NA", "primary", "metastasis")
grouplist <- factor(grouplist, levels = c("metastasis", "primary"))

save(exp, ID, ph, file = "GSE41568_Raw_Data.Rdata")
save(grouplist, file = "group_metastasis_primary.Rdata")

# ---- GSE71222 ----
setwd("path_to/GSE71222")
eSet <- getGEO(filename = "GSE71222_series_matrix.txt.gz", getGPL = FALSE)
exp <- exprs(eSet)
exp <- normalize_expression(exp)

ph <- as.data.frame(pData(eSet))
ID <- as.data.frame(read_excel("GPL570.xlsx", sheet = 1, col_names = TRUE))
exp <- map_gene_symbols(exp, ID, "Gene Symbol")
ph <- ph[ph$`tissue:ch1` %in% c("cancer, homogenized ", "cancer, LCM"), ]
tmp <- filter_samples(ph, exp, "rownames(ph)")
ph <- tmp$ph; exp <- tmp$exp
rownames(exp) <- make.unique(sub(" /// .*", "", rownames(exp)))
grouplist <- ifelse(ph$`metastatic tumor site:ch1` == "NA", "primary", "metastasis")
grouplist <- factor(grouplist, levels = c("metastasis", "primary"))

save(exp, ID, ph, file = "GSE71222_Raw_Data.Rdata")
save(grouplist, file = "group_metastasis_primary.Rdata")

# ---- GSE81986 ----
setwd("path_to/GSE81986")
eSet <- getGEO(filename = "GSE81986-GPL570_series_matrix (1).txt.gz", getGPL = FALSE)
exp <- exprs(eSet)
exp <- normalize_expression(exp)

ph <- as.data.frame(pData(eSet))
ID <- as.data.frame(read_excel("GPL570.xlsx", sheet = 1, col_names = TRUE))
exp <- map_gene_symbols(exp, ID, "Gene Symbol")
tmp <- filter_samples(ph, exp, "rownames(ph)")
ph <- tmp$ph; exp <- tmp$exp
rownames(exp) <- make.unique(sub(" /// .*", "", rownames(exp)))
grouplist <- ifelse(ph$`metastatic tumor site:ch1` == "NA", "primary", "metastasis")
grouplist <- factor(grouplist, levels = c("metastasis", "primary"))

save(exp, ID, ph, file = "GSE81986_Adjusted_Data_Gene_Symbol.Rdata")
save(grouplist, file = "group_metastasis_primary.Rdata")

# ============================================================================
# End of Preprocessing Pipeline
# ============================================================================
