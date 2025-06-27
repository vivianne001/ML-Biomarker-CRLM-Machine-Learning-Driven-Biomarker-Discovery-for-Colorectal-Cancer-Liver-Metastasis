# Regularization and Feature Subset Selection for Machine Learning Inputs

# Load required packages
library(mlr3)
library(mlr3extralearners)
library(mlr3verse)
library(dplyr)
library(tidyverse)

# View available classification learners
learner_list <- as.data.table(mlr_learners)
classification_learners <- learner_list[task_type == "classif"]
save(learner_list, file = "mlr3_classification_learners_list.Rdata")

# 1. Organize data frames for analysis

# Prepare training set - filter and transpose as needed
k <- rownames(trainset)
trainset <- exp[, k]
trainset <- trainset[dd, ]
trainset <- na.omit(trainset)

# Prepare validation set - filter and transpose as needed
k <- rownames(validationset)
validationset <- exp[, k]
validationset <- validationset[dd, ]
validationset <- na.omit(validationset)
ph_valid <- ph[k, ]
validationset <- t(validationset)
validationset <- as.data.frame(validationset)
validationset$status <- ph_valid$group
validationset$status <- recode(validationset$status, "primary" = 0, "metastasis" = 1)
validationset$status <- as.factor(validationset$status)

# Prepare test set - filter and transpose as needed
k <- rownames(testset)
testset <- exp[, k]
testset <- testset[dd, ]
testset <- na.omit(testset)
ph_test <- ph[k, ]
testset <- t(testset)
testset <- as.data.frame(testset)
testset$status <- ph_test$group
testset$status <- recode(testset$status, "primary" = 0, "metastasis" = 1)
testset$status <- as.factor(testset$status)

# Final preparation of training set
trainset <- t(trainset)
trainset <- as.data.frame(trainset)
ph_train <- ph[rownames(trainset), ]

index <- intersect(rownames(ph_train), rownames(trainset))
trainset <- trainset[index, ]
ph_train <- ph_train[index, ]
trainset$status <- ph_train$group

# Convert status column to numeric or logical
table(trainset$status)
class(trainset$status)

# Convert status labels to 0 and 1
trainset$status <- recode(trainset$status, "primary" = 0, "metastasis" = 1)
table(trainset$status)

mpe <- trainset
mpe$status <- as.factor(mpe$status)

# Fix column names to comply with R variable naming conventions
invalid_colnames <- colnames(mpe)[!make.names(colnames(mpe), unique = TRUE) %in% colnames(mpe)]
print(invalid_colnames)
new_colnames <- make.names(colnames(mpe), unique = TRUE)
colnames(mpe) <- new_colnames
colnames(mpe) <- gsub("^(.*?)\\..*$", "\\1", colnames(mpe))

# Select relevant columns (make sure 'cols' is defined elsewhere)
mpe <- mpe[, cols]
# Set up classification task
task <- as_task_classif(mpe, target = "status")

# 2. Feature selection using wrapper method: forward stepwise selection

library(mlr3fselect)
set.seed(123)

learner <- lrn("classif.cv_glmnet")
resampling <- rsmp("cv", folds = 10)
resampling$instantiate(task)

learner$train(task)
filter <- flt("selected_features", learner = learner)
filter$calculate(task)
flt <- as.data.table(filter)

save(flt, file = "feature_filter_set_seed.Rdata")

# Feature Selection forming feature subset
flt = flt[flt$score = 1,]