# ============================================================================
# Combine and Split 6 GEO Datasets into Training, Validation, and Test Sets
# ============================================================================

rm(list = ls())

# ===============================
# Load Preprocessed Data
# ===============================

# Please set the correct paths to your .Rdata files
load("GSE81986_Raw_Data.Rdata")   # exp, ID, ph
load("group_metastasis_primary_GSE81986.Rdata") # grouplist
exp1 <- exp; ph1 <- ph; group1 <- grouplist; dataset1 <- "GSE81986"

load("GSE27854_Raw_Data.Rdata")
load("group_metastasis_primary_GSE27854.Rdata")
exp2 <- exp; ph2 <- ph; group2 <- grouplist; dataset2 <- "GSE27854"

load("GSE71222_Raw_Data.Rdata")
load("group_metastasis_primary_GSE71222.Rdata")
exp3 <- exp; ph3 <- ph; group3 <- grouplist; dataset3 <- "GSE71222"

load("GSE41568_Raw_Data.Rdata")
load("group_metastasis_primary_GSE41568.Rdata")
exp4 <- exp; ph4 <- ph; group4 <- grouplist; dataset4 <- "GSE41568"

load("GSE51244_Raw_Data.Rdata")
load("group_metastasis_primary_GSE51244.Rdata")
exp5 <- exp; ph5 <- ph; group5 <- grouplist; dataset5 <- "GSE51244"

load("GSE21510_Raw_Data.Rdata")
load("group_metastasis_primary_GSE21510.Rdata")
exp6 <- exp; ph6 <- ph; group6 <- grouplist; dataset6 <- "GSE21510"

load("GSE18105_Raw_Data.Rdata")
load("group_metastasis_primary_GSE18105.Rdata")
exp7 <- exp; ph7 <- ph; group7 <- grouplist; dataset7 <- "GSE18105"

# ===============================
# Prepare Expression and Phenotype Data
# ===============================

# Add sample IDs and dataset labels
exp1$ID <- rownames(exp1); ph1$group <- group1; ph1$dataset <- dataset1
exp2$ID <- rownames(exp2); ph2$group <- group2; ph2$dataset <- dataset2
exp3$ID <- rownames(exp3); ph3$group <- group3; ph3$dataset <- dataset3
exp4$ID <- rownames(exp4); ph4$group <- group4; ph4$dataset <- dataset4
exp5$ID <- rownames(exp5); ph5$group <- group5; ph5$dataset <- dataset5
exp6$ID <- rownames(exp6); ph6$group <- group6; ph6$dataset <- dataset6
exp7$ID <- rownames(exp7); ph7$group <- group7; ph7$dataset <- dataset7

# ===============================
# Merge Expression Matrices by Gene ID
# ===============================

exp_merged <- Reduce(function(x, y) merge(x, y, by = "ID"), 
                     list(exp1, exp2, exp3, exp4, exp5, exp6, exp7))
rownames(exp_merged) <- exp_merged$ID
exp_merged <- exp_merged[ , -1] # Remove ID column

# ===============================
# Combine Phenotype Data
# ===============================

ph_combined <- rbind(
  ph1, ph2, ph3, ph4, ph5, ph6, ph7
)

# ===============================
# Assign Datasets to Train, Validation, and Test Sets
# ===============================

# Define which datasets go to which set
train_datasets <- c("GSE18105", "GSE81986", "GSE41568", "GSE71222")
validation_datasets <- c("GSE51244")
test_datasets <- c("GSE21510")

# Training set
train_samples <- ph_combined[ph_combined$dataset %in% train_datasets, "geo_accession"]
trainset <- t(exp_merged[, train_samples])
trainset <- as.data.frame(trainset)
trainset$status <- ph_combined[match(rownames(trainset), ph_combined$geo_accession), "group"]
trainset$status <- factor(ifelse(trainset$status == "primary", 0, 1), levels = c(0, 1))

# Validation set
validation_samples <- ph_combined[ph_combined$dataset %in% validation_datasets, "geo_accession"]
validationset <- t(exp_merged[, validation_samples])
validationset <- as.data.frame(validationset)
validationset$status <- ph_combined[match(rownames(validationset), ph_combined$geo_accession), "group"]
validationset$status <- factor(ifelse(validationset$status == "primary", 0, 1), levels = c(0, 1))

# Test set
test_samples <- ph_combined[ph_combined$dataset %in% test_datasets, "geo_accession"]
testset <- t(exp_merged[, test_samples])
testset <- as.data.frame(testset)
testset$status <- ph_combined[match(rownames(testset), ph_combined$geo_accession), "group"]
testset$status <- factor(ifelse(testset$status == "primary", 0, 1), levels = c(0, 1))

# ===============================
# Save Results
# ===============================
save(trainset, validationset, testset, file = "train_validation_test_sets.Rdata")

# ============================================================================
# End of Script
# ============================================================================
