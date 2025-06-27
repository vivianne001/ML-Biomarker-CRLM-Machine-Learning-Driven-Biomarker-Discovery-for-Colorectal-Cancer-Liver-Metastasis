<<<<<<< HEAD
# Comprehensive Classification Evaluation for Multiple Machine Learning Models: ACC, AUC, Confusion Matrix, and ROC Curve Visualization

rm(list = ls())

# Load necessary libraries
library(mlr3)
library(mlr3learners)
library(mlr3extralearners)
library(mlr3filters)
library(mlr3verse)
library(mlr3viz)
library(ggplot2)
library(pROC)

# Ensure 'mpe', 'ph', 'trainset', 'testset', 'validationset' are loaded and properly formatted

# Standardize column names
colnames(mpe) <- gsub("^(.*?)\\..*$", "\\1", colnames(mpe))
mpe$status <- as.factor(mpe$status)
ph$group <- ifelse(ph$group == "metastasis", 1, 0)
ph$group <- as.factor(ph$group)

# Create mlr3 tasks for each dataset
task <- as_task_classif(trainset, target = "status")
test_task <- as_task_classif(testset, target = "status")
validation_task <- as_task_classif(validationset, target = "status")

# Example: Evaluate kknn classifier (add more learners for multiple models)
learner <- lrn("classif.kknn", predict_type = "prob", k = 7)
resampling <- rsmp("cv", folds = 10)
resampling$instantiate(task)

# Other models for this classification task
# classif.debug
# classif.earth
# classif.featureless
# classif.C50
# classif.gam
# classif.gamboost
# classif.gausspr
# classif.gbm
# classif.glmboost
# classif.glmnet
# classif.knn
# classif.lda
# classif.lightgbm
# classif.log_reg
# classif.naivebayes
# classif.rpart
# classif.svm
# classif.xgboost
# classif.multinom
# classif.ranger

# Train the model
learner$train(task)

# Evaluate on validation set
validation_prediction <- learner$predict(validation_task)
print(validation_prediction)
confusion_validation <- validation_prediction$confusion
save(confusion_validation, file = "confusion_matrix_validation.Rdata")

# Evaluate on test set
test_prediction <- learner$predict(test_task)
confusion_test <- test_prediction$confusion
save(confusion_test, file = "confusion_matrix_test.Rdata")

# Compute accuracy and AUC for validation set
acc_validation <- validation_prediction$score(msr("classif.acc"))
auc_validation <- validation_prediction$score(msr("classif.auc"))
save(acc_validation, auc_validation, file = "acc_auc_validation.Rdata")

# Compute accuracy and AUC for test set
acc_test <- test_prediction$score(msr("classif.acc"))
auc_test <- test_prediction$score(msr("classif.auc"))
save(acc_test, auc_test, file = "acc_auc_test.Rdata")

# Visualize accuracy results on validation set
autoplot(validation_prediction, measure = msr("classif.acc")) +
  scale_fill_manual(values = c("#EF8A43", "#4865A9"), name = "Group", labels = c("primary CRC", "CRLM")) +
  labs(title = "Model Accuracy on Validation Dataset", y = "Sample Count") +
  scale_x_discrete(labels = c("Actual Grouping", "Predictive Grouping")) +
  theme(axis.line = element_line(color = "black"))

# ROC curve for validation set
auc_val <- validation_prediction$score(msr("classif.auc"))
autoplot(validation_prediction, type = "roc") +
  geom_text(aes(label = paste("Validation Dataset AUC =", round(auc_val, 3))),
            x = .6, y = .1, hjust = 0, vjust = .5, size = 4, color = "#345D82") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggtitle("ROC Curve of Validation Dataset") +
  theme(axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

# ROC curve for test set
auc_test_val <- test_prediction$score(msr("classif.auc"))
autoplot(test_prediction, type = "roc") +
  geom_text(aes(label = paste("Test Dataset AUC =", round(auc_test_val, 3))),
            x = .6, y = .1, hjust = 0, vjust = .5, size = 4.3, color = "#345D82") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggtitle("ROC Curve of Test Dataset") +
  theme(axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

# Additional ROC curve using pROC for validation set
probabilities_val <- validation_prediction$prob[, 2] # Assuming class "1" is positive
actual_labels_val <- validation_prediction$truth
roc_data_val <- roc(actual_labels_val, probabilities_val, levels = c("0", "1"))
roc_df_val <- data.frame(
  specificity = roc_data_val$specificities,
  sensitivity = roc_data_val$sensitivities
)
auc_value_val <- auc(roc_data_val)
auc_ci_val <- ci.auc(roc_data_val)
auc_text_val <- sprintf("AUC = %.4f (95%% CI: %.2f-%.2f)", auc_value_val, auc_ci_val[1], auc_ci_val[3])
ggplot(roc_df_val, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line() +
  geom_abline(linetype = "dashed", color = "grey") +
  scale_x_continuous(limits = c(0, 1), name = "1-Specificity") +
  scale_y_continuous(limits = c(0, 1), name = "Sensitivity") +
  ggtitle("ROC Curve of Validation Dataset") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black")) +
  annotate("text", x = 0.6, y = 0.15, label = "mCRC Liver Metastasis & Primary", size = 4, hjust = 0) +
  annotate("text", x = 0.6, y = 0.1, label = auc_text_val, size = 4, hjust = 0)

# You can repeat similar steps for other models or for more advanced visualization.
=======
# Comprehensive Classification Evaluation for Multiple Machine Learning Models: ACC, AUC, Confusion Matrix, and ROC Curve Visualization

rm(list = ls())

# Load necessary libraries
library(mlr3)
library(mlr3learners)
library(mlr3extralearners)
library(mlr3filters)
library(mlr3verse)
library(mlr3viz)
library(ggplot2)
library(pROC)

# Ensure 'mpe', 'ph', 'trainset', 'testset', 'validationset' are loaded and properly formatted

# Standardize column names
colnames(mpe) <- gsub("^(.*?)\\..*$", "\\1", colnames(mpe))
mpe$status <- as.factor(mpe$status)
ph$group <- ifelse(ph$group == "metastasis", 1, 0)
ph$group <- as.factor(ph$group)

# Create mlr3 tasks for each dataset
task <- as_task_classif(trainset, target = "status")
test_task <- as_task_classif(testset, target = "status")
validation_task <- as_task_classif(validationset, target = "status")

# Example: Evaluate kknn classifier (add more learners for multiple models)
learner <- lrn("classif.kknn", predict_type = "prob", k = 7)
resampling <- rsmp("cv", folds = 10)
resampling$instantiate(task)

# Other models for this classification task
# classif.debug
# classif.earth
# classif.featureless
# classif.C50
# classif.gam
# classif.gamboost
# classif.gausspr
# classif.gbm
# classif.glmboost
# classif.glmnet
# classif.knn
# classif.lda
# classif.lightgbm
# classif.log_reg
# classif.naivebayes
# classif.rpart
# classif.svm
# classif.xgboost
# classif.multinom
# classif.ranger

# Train the model
learner$train(task)

# Evaluate on validation set
validation_prediction <- learner$predict(validation_task)
print(validation_prediction)
confusion_validation <- validation_prediction$confusion
save(confusion_validation, file = "confusion_matrix_validation.Rdata")

# Evaluate on test set
test_prediction <- learner$predict(test_task)
confusion_test <- test_prediction$confusion
save(confusion_test, file = "confusion_matrix_test.Rdata")

# Compute accuracy and AUC for validation set
acc_validation <- validation_prediction$score(msr("classif.acc"))
auc_validation <- validation_prediction$score(msr("classif.auc"))
save(acc_validation, auc_validation, file = "acc_auc_validation.Rdata")

# Compute accuracy and AUC for test set
acc_test <- test_prediction$score(msr("classif.acc"))
auc_test <- test_prediction$score(msr("classif.auc"))
save(acc_test, auc_test, file = "acc_auc_test.Rdata")

# Visualize accuracy results on validation set
autoplot(validation_prediction, measure = msr("classif.acc")) +
  scale_fill_manual(values = c("#EF8A43", "#4865A9"), name = "Group", labels = c("primary CRC", "CRLM")) +
  labs(title = "Model Accuracy on Validation Dataset", y = "Sample Count") +
  scale_x_discrete(labels = c("Actual Grouping", "Predictive Grouping")) +
  theme(axis.line = element_line(color = "black"))

# ROC curve for validation set
auc_val <- validation_prediction$score(msr("classif.auc"))
autoplot(validation_prediction, type = "roc") +
  geom_text(aes(label = paste("Validation Dataset AUC =", round(auc_val, 3))),
            x = .6, y = .1, hjust = 0, vjust = .5, size = 4, color = "#345D82") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggtitle("ROC Curve of Validation Dataset") +
  theme(axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

# ROC curve for test set
auc_test_val <- test_prediction$score(msr("classif.auc"))
autoplot(test_prediction, type = "roc") +
  geom_text(aes(label = paste("Test Dataset AUC =", round(auc_test_val, 3))),
            x = .6, y = .1, hjust = 0, vjust = .5, size = 4.3, color = "#345D82") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggtitle("ROC Curve of Test Dataset") +
  theme(axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

# Additional ROC curve using pROC for validation set
probabilities_val <- validation_prediction$prob[, 2] # Assuming class "1" is positive
actual_labels_val <- validation_prediction$truth
roc_data_val <- roc(actual_labels_val, probabilities_val, levels = c("0", "1"))
roc_df_val <- data.frame(
  specificity = roc_data_val$specificities,
  sensitivity = roc_data_val$sensitivities
)
auc_value_val <- auc(roc_data_val)
auc_ci_val <- ci.auc(roc_data_val)
auc_text_val <- sprintf("AUC = %.4f (95%% CI: %.2f-%.2f)", auc_value_val, auc_ci_val[1], auc_ci_val[3])
ggplot(roc_df_val, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line() +
  geom_abline(linetype = "dashed", color = "grey") +
  scale_x_continuous(limits = c(0, 1), name = "1-Specificity") +
  scale_y_continuous(limits = c(0, 1), name = "Sensitivity") +
  ggtitle("ROC Curve of Validation Dataset") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black")) +
  annotate("text", x = 0.6, y = 0.15, label = "mCRC Liver Metastasis & Primary", size = 4, hjust = 0) +
  annotate("text", x = 0.6, y = 0.1, label = auc_text_val, size = 4, hjust = 0)

# You can repeat similar steps for other models or for more advanced visualization.
>>>>>>> 027b93a2ffda7f6e9e9524ea22593a1cf937035a
# Multiple models can be compared by looping over different learners and saving all performance metrics.