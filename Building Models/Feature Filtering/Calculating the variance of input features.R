# Calculating the variance of input features.

# Load necessary libraries
library(mlr3verse)
library(dplyr)
library(data.table)
library(FSelectorRcpp)

# Ensure the necessary package is installed
install.packages('FSelectorRcpp')

# Prepare the training set with functional genes
index <- intersect(colnames(trainset), function_genes$geneSymbols)
index <- c(index, "status")
trainset <- trainset[, index]

# Create a classification task
task <- as_task_classif(trainset, target = "status")

# Apply variance filter to select features
filter <- flt("variance")
filter$calculate(task)
filt <- as.data.table(filter)

# Save the filtering result
save(filt, file = "filter-variance.Rdata")

# Interpretation of Scores:
# High scores indicate features with large variance and greater contribution to model performance.
# Low scores indicate features with small variance, and they may be redundant or contribute less to the model.