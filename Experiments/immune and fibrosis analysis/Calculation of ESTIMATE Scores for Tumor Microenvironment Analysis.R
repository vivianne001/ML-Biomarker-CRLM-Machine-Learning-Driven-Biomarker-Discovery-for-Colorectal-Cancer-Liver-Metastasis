################################################################################
# Calculation of ESTIMATE Scores for Tumor Microenvironment Analysis
################################################################################

# Section 1: Setup and Environment Preparation

rm(list = ls())

library(IOBR)
library(EPIC)
library(estimate)
library(tidyverse)
library(tidyHeatmap)
library(maftools)
library(ggpubr)
library(ggplot2)
library(reshape)

# Section 2: Data Loading and Preprocessing

# Example: load("your_data.Rdata")
# Select study cohort (High & Low drug sensitive groups)
k <- new_data$id
exp <- exp[, k]
clin <- clin[k, ]

head(exp)
range(exp)

# Section 3: Counts to TPM Conversion

options(stringsAsFactors = FALSE)
library(dplyr)

# Prepare gene annotation for TPM conversion
ann <- exp
ann$Gene <- str_split(string = rownames(ann), pattern = " ", simplify = TRUE)[, 1]
input <- read.csv(file = "All_hg19gene_len.csv")
merge <- left_join(ann, input, by = "Gene")
merge <- na.omit(merge)
merge[1:10, 120:124]

# Calculate TPM
kb <- merge$Length / 1000
countdata <- merge[, 1:119]
rpk <- countdata / kb
tpm <- t(t(rpk) / colSums(rpk) * 1e6)

range(exp)
range(tpm)

tpm <- as.data.frame(tpm)
rownames(tpm) <- merge$Gene
tpm <- tpm[!duplicated(rownames(tpm)), ]
tpm <- log2(tpm + 1)
range(tpm)
save(tpm, file = "exp of tpm.Rdata")

# Section 4: ESTIMATE Score Calculation

estimate.normalized <- deconvo_tme(eset = tpm, method = "estimate")
TCGA_TME.results <- estimate.normalized

# save(TCGA_TME.results, file = "TCGA_LIHC_TME.results.Rdata")
save(estimate.normalized, file = "estimate.tpm.Rdata")

# Section 5: Group Information Preparation

group_list <- new_data$new_column
group_list <- as.factor(group_list)
table(group_list)