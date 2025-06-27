# ============================================================================
# Merge Four GEO Datasets into Training Set and Visualize PCA Before/After Batch Effect Removal
# ============================================================================

rm(list = ls())

# Load required libraries
library(limma)
library(sva)
library(ggplot2)

# Load preprocessed data for each dataset
load("GSE18105_Raw_Data.Rdata")        # exp, ID, ph
load("group_metastasis_primary_GSE18105.Rdata") # grouplist
exp1 <- exp; ph1 <- ph; group1 <- grouplist; dataset1 <- "GSE18105"

load("GSE81986_Raw_Data.Rdata")
load("group_metastasis_primary_GSE81986.Rdata")
exp2 <- exp; ph2 <- ph; group2 <- grouplist; dataset2 <- "GSE81986"

load("GSE41568_Raw_Data.Rdata")
load("group_metastasis_primary_GSE41568.Rdata")
exp3 <- exp; ph3 <- ph; group3 <- grouplist; dataset3 <- "GSE41568"

load("GSE71222_Raw_Data.Rdata")
load("group_metastasis_primary_GSE71222.Rdata")
exp4 <- exp; ph4 <- ph; group4 <- grouplist; dataset4 <- "GSE71222"

# Add gene ID columns and dataset labels
exp1$ID <- rownames(exp1); ph1$group <- group1; ph1$dataset <- dataset1
exp2$ID <- rownames(exp2); ph2$group <- group2; ph2$dataset <- dataset2
exp3$ID <- rownames(exp3); ph3$group <- group3; ph3$dataset <- dataset3
exp4$ID <- rownames(exp4); ph4$group <- group4; ph4$dataset <- dataset4

# Merge expression matrices by gene ID
exp_merged <- Reduce(function(x, y) merge(x, y, by = "ID"), 
                     list(exp1, exp2, exp3, exp4))
rownames(exp_merged) <- exp_merged$ID
exp_merged <- exp_merged[ , -1]

# Merge phenotype data
ph_merged <- rbind(ph1, ph2, ph3, ph4)

# Extract training set sample names
train_samples <- ph_merged$geo_accession

# Create training set
trainset <- t(exp_merged[, train_samples])
trainset <- as.data.frame(trainset)
trainset$status <- ph_merged[match(rownames(trainset), ph_merged$geo_accession), "group"]
trainset$status <- factor(ifelse(trainset$status == "primary", 0, 1), levels = c(0, 1))

# ===============================
# PCA Visualization Before Batch Effect Removal
# ===============================

# Prepare data for PCA (remove status column)
trainset_numeric <- trainset[, !(colnames(trainset) %in% "status")]
batch_info <- ph_merged[match(rownames(trainset), ph_merged$geo_accession), "dataset"]

# PCA before batch effect removal
pca_before <- prcomp(trainset_numeric, scale. = TRUE)
pca_data_before <- data.frame(pca_before$x, Batch = batch_info, Status = trainset$status)

# Plot PCA before batch effect removal
ggplot(pca_data_before, aes(x = PC1, y = PC2, color = Batch, shape = Status)) +
  geom_point(size = 2) +
  ggtitle("PCA Before Batch Effect Removal") +
  theme_minimal()

# ===============================
# Batch Effect Removal
# ===============================

# Remove batch effect using limma
design <- model.matrix(~ trainset$status)
trainset_corrected <- t(removeBatchEffect(t(trainset_numeric), batch = batch_info, design = design))

# Alternatively, use ComBat (uncomment if you prefer ComBat)
# trainset_corrected <- ComBat(dat = t(trainset_numeric), batch = batch_info)
# trainset_corrected <- t(trainset_corrected)

# ===============================
# PCA Visualization After Batch Effect Removal
# ===============================

pca_after <- prcomp(trainset_corrected, scale. = TRUE)
pca_data_after <- data.frame(pca_after$x, Batch = batch_info, Status = trainset$status)

ggplot(pca_data_after, aes(x = PC1, y = PC2, color = Batch, shape = Status)) +
  geom_point(size = 2) +
  ggtitle("PCA After Batch Effect Removal") +
  theme_minimal()

# ===============================
# Save Results
# ===============================

save(trainset_corrected, file = "removebatcheffects.trainset.Rdata")

# ============================================================================
# End of Script
# ============================================================================
