<<<<<<< HEAD
# Heatmap Visualization of AUC Values for 29 Candidate Biomarkers

# Load AUC values for validation set
load("feature_auc_values-testset validation-permutation-trainset 29.Rdata")
validation_auc <- as.data.frame(feature_auc_values)
validation_auc <- as.data.frame(t(validation_auc))
colnames(validation_auc) <- c("Validation_set")
validation_auc$genes <- rownames(validation_auc)

# Load AUC values for training set
load("feature_auc_values-permutation-trainset 29 + CRC BIOMARKER Features.Rdata")
trainset_auc <- as.data.frame(feature_auc_values)
trainset_auc <- as.data.frame(t(trainset_auc))
colnames(trainset_auc) <- c("Training_set")

# Load selected feature names
load("特征过滤-set.seed.Rdata")
flt <- flt[flt$score == 1, ]
genes <- flt$feature

trainset_auc$genes <- rownames(trainset_auc)
trainset_auc <- trainset_auc[genes, ]

# Merge training and validation AUC matrices
auc <- merge(trainset_auc, validation_auc, by = "genes")

# Order features by training set AUC
library(dplyr)
auc <- arrange(auc, desc(Training_set))
rownames(auc) <- auc$genes
auc <- auc[, -1]

# Set new column names
colnames(auc) <- c("Training_set", "Validation_set")

# Create a heatmap for AUC values
library(ComplexHeatmap)
library(circlize)

auc_matrix <- as.matrix(auc)
col_fun <- colorRamp2(c(0.73, 0.83, 0.93), c("#4878b5", "#fefec1", "#d83128"))

Heatmap(
  auc_matrix,
  name = "AUC",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  column_gap = unit(2, "mm"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 12),
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.4f", auc_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
  }
=======
# Heatmap Visualization of AUC Values for 29 Candidate Biomarkers

# Load AUC values for validation set
load("feature_auc_values-testset validation-permutation-trainset 29.Rdata")
validation_auc <- as.data.frame(feature_auc_values)
validation_auc <- as.data.frame(t(validation_auc))
colnames(validation_auc) <- c("Validation_set")
validation_auc$genes <- rownames(validation_auc)

# Load AUC values for training set
load("feature_auc_values-permutation-trainset 29 + CRC BIOMARKER Features.Rdata")
trainset_auc <- as.data.frame(feature_auc_values)
trainset_auc <- as.data.frame(t(trainset_auc))
colnames(trainset_auc) <- c("Training_set")

# Load selected feature names
load("特征过滤-set.seed.Rdata")
flt <- flt[flt$score == 1, ]
genes <- flt$feature

trainset_auc$genes <- rownames(trainset_auc)
trainset_auc <- trainset_auc[genes, ]

# Merge training and validation AUC matrices
auc <- merge(trainset_auc, validation_auc, by = "genes")

# Order features by training set AUC
library(dplyr)
auc <- arrange(auc, desc(Training_set))
rownames(auc) <- auc$genes
auc <- auc[, -1]

# Set new column names
colnames(auc) <- c("Training_set", "Validation_set")

# Create a heatmap for AUC values
library(ComplexHeatmap)
library(circlize)

auc_matrix <- as.matrix(auc)
col_fun <- colorRamp2(c(0.73, 0.83, 0.93), c("#4878b5", "#fefec1", "#d83128"))

Heatmap(
  auc_matrix,
  name = "AUC",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  column_gap = unit(2, "mm"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 12),
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.4f", auc_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
  }
>>>>>>> 027b93a2ffda7f6e9e9524ea22593a1cf937035a
)