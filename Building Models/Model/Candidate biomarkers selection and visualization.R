<<<<<<< HEAD
######### AUC Calculation and Visualization for 29 Candidate Biomarkers #########

rm(list = ls())

library(mlr3)
library(mlr3learners)
library(mlr3measures)
library(pROC)
library(ggplot2)
library(ROCR)
library(ComplexHeatmap)
library(circlize)

# Load preprocessed datasets and filtered features
load("removebatcheffects.trainset-all genes plus status.Rdata")
flt = flt[flt$score == 1,]
feature = flt$feature

# Subset to the required features and status
trainset = trainset[,c("CXCL12","CXCR4","EGFR","PTGS2","ERBB2","ERBB3","CDX2","FGF19",feature,"CDH17","YBX1","METTL14","CTSK","TERT","KRT18","KRT15","KRT7","GPA33","SATB2","status")]

# Initialize lists to store AUC values and prediction probabilities
feature_predictions <- list()
feature_auc_values <- list()

# Calculate AUC for each biomarker
for (feature_name in colnames(trainset)[-ncol(trainset)]) {
  single_feature_trainset <- trainset[, c(feature_name, "status"), drop = FALSE]
  single_feature_task <- as_task_classif(single_feature_trainset, target = "status")
  single_feature_test = as_task_classif(testset, target = "status")
  learner <- lrn("classif.kknn", predict_type = "prob")
  learner$train(single_feature_task)
  prediction <- learner$predict(single_feature_test)
  auc_measure <- msr("classif.auc")
  auc_value <- auc_measure$score(prediction)
  feature_auc_values[[feature_name]] <- auc_value
  feature_predictions[[feature_name]] <- prediction$prob[,2]
}

# Convert AUC values to data frame
auc_values = as.data.frame(feature_auc_values)
auc_value = t(auc_values)
auc_value = as.data.frame(auc_value)

# Save results
save(feature_auc_values, feature_predictions, file = "feature_auc_values-permutation-trainset 29 + CRC BIOMARKER Features.Rdata")
save(feature_auc_values, feature_predictions, file = "feature_auc_values-testset validation-permutation-trainset 29.Rdata")

# Visualize ROC curves for a subset of genes
genes <- c("ACMSD","TERT","KRT18","KRT15","KRT7","CDH17","GPA33","SATB2","CDX2")
colors <- c("ACMSD" = "#B42B22", "CDX2" = "#EC3232", "SATB2" = "#F5CBBF", "GPA33" = "#0787C3","CDH17" = "#D5EAED","KRT7" = "#996E2E","KRT15" = "#F6944B","KRT18" = "#F4DBB2","TERT" = "#315A89")

roc_df <- data.frame()
auc_infos <- c()

trainset1$status = ifelse(trainset1$status == "metastasis",1,0)
trainset1$status = as.factor(trainset1$status)

additional_texts <- c("TERT" = "Bertorelle,R et al.2013","KRT18" = "Zhang,Jingfeng et al.2019","KRT15" = "Rao,X et al.2020",
                      "KRT7" = "Czapiewski,Piotr et al.2016","CDH17" = "Brandler,Tamar C et al.2015","GPA33" = "Wong,Newton A C S et al.2017","SATB2" = "Neri,Giuseppe et al.2020",
                      "CDX2" = "Toth,Csaba et al.2018")

for(gene in genes) {
  roc_data <- roc(response = trainset1$status, predictor = feature_predictions[[gene]], levels=c("1","0"))
  temp_df <- data.frame(
    specificity = roc_data$specificities,
    sensitivity = roc_data$sensitivities,
    gene = gene
  )
  roc_df <- rbind(roc_df, temp_df)
  auc_value <- auc(roc_data)
  ci <- ci.auc(roc_data)
  specific_text <- additional_texts[gene]
  auc_info <- sprintf("%s: AUC = %.4f (95%% CI: %.4f-%.4f) %s", gene, auc_value, ci[1], ci[3], specific_text)
  auc_infos <- c(auc_infos, auc_info)
}

p <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, color = gene)) +
  geom_line() +
  geom_abline(linetype = "dashed", color = "grey") +
  scale_x_continuous(limits = c(0, 1), name = "1-Specificity") +
  scale_y_continuous(limits = c(0, 1), name = "Sensitivity") +
  ggtitle("ROC curve for multiple genes") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black")) +
  scale_color_manual(values = colors)

for(i in 1:length(auc_infos)) {
  p <- p + annotate("text", x = 0.55, y = 0.1 + (i-1)*0.05, label = auc_infos[i], size = 3, hjust = 0, color = "black")
}

print(p)
=======
######### AUC Calculation and Visualization for 29 Candidate Biomarkers #########

rm(list = ls())

library(mlr3)
library(mlr3learners)
library(mlr3measures)
library(pROC)
library(ggplot2)
library(ROCR)
library(ComplexHeatmap)
library(circlize)

# Load preprocessed datasets and filtered features
load("removebatcheffects.trainset-all genes plus status.Rdata")
flt = flt[flt$score == 1,]
feature = flt$feature

# Subset to the required features and status
trainset = trainset[,c("CXCL12","CXCR4","EGFR","PTGS2","ERBB2","ERBB3","CDX2","FGF19",feature,"CDH17","YBX1","METTL14","CTSK","TERT","KRT18","KRT15","KRT7","GPA33","SATB2","status")]

# Initialize lists to store AUC values and prediction probabilities
feature_predictions <- list()
feature_auc_values <- list()

# Calculate AUC for each biomarker
for (feature_name in colnames(trainset)[-ncol(trainset)]) {
  single_feature_trainset <- trainset[, c(feature_name, "status"), drop = FALSE]
  single_feature_task <- as_task_classif(single_feature_trainset, target = "status")
  single_feature_test = as_task_classif(testset, target = "status")
  learner <- lrn("classif.kknn", predict_type = "prob")
  learner$train(single_feature_task)
  prediction <- learner$predict(single_feature_test)
  auc_measure <- msr("classif.auc")
  auc_value <- auc_measure$score(prediction)
  feature_auc_values[[feature_name]] <- auc_value
  feature_predictions[[feature_name]] <- prediction$prob[,2]
}

# Convert AUC values to data frame
auc_values = as.data.frame(feature_auc_values)
auc_value = t(auc_values)
auc_value = as.data.frame(auc_value)

# Save results
save(feature_auc_values, feature_predictions, file = "feature_auc_values-permutation-trainset 29 + CRC BIOMARKER Features.Rdata")
save(feature_auc_values, feature_predictions, file = "feature_auc_values-testset validation-permutation-trainset 29.Rdata")

# Visualize ROC curves for a subset of genes
genes <- c("ACMSD","TERT","KRT18","KRT15","KRT7","CDH17","GPA33","SATB2","CDX2")
colors <- c("ACMSD" = "#B42B22", "CDX2" = "#EC3232", "SATB2" = "#F5CBBF", "GPA33" = "#0787C3","CDH17" = "#D5EAED","KRT7" = "#996E2E","KRT15" = "#F6944B","KRT18" = "#F4DBB2","TERT" = "#315A89")

roc_df <- data.frame()
auc_infos <- c()

trainset1$status = ifelse(trainset1$status == "metastasis",1,0)
trainset1$status = as.factor(trainset1$status)

additional_texts <- c("TERT" = "Bertorelle,R et al.2013","KRT18" = "Zhang,Jingfeng et al.2019","KRT15" = "Rao,X et al.2020",
                      "KRT7" = "Czapiewski,Piotr et al.2016","CDH17" = "Brandler,Tamar C et al.2015","GPA33" = "Wong,Newton A C S et al.2017","SATB2" = "Neri,Giuseppe et al.2020",
                      "CDX2" = "Toth,Csaba et al.2018")

for(gene in genes) {
  roc_data <- roc(response = trainset1$status, predictor = feature_predictions[[gene]], levels=c("1","0"))
  temp_df <- data.frame(
    specificity = roc_data$specificities,
    sensitivity = roc_data$sensitivities,
    gene = gene
  )
  roc_df <- rbind(roc_df, temp_df)
  auc_value <- auc(roc_data)
  ci <- ci.auc(roc_data)
  specific_text <- additional_texts[gene]
  auc_info <- sprintf("%s: AUC = %.4f (95%% CI: %.4f-%.4f) %s", gene, auc_value, ci[1], ci[3], specific_text)
  auc_infos <- c(auc_infos, auc_info)
}

p <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, color = gene)) +
  geom_line() +
  geom_abline(linetype = "dashed", color = "grey") +
  scale_x_continuous(limits = c(0, 1), name = "1-Specificity") +
  scale_y_continuous(limits = c(0, 1), name = "Sensitivity") +
  ggtitle("ROC curve for multiple genes") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black")) +
  scale_color_manual(values = colors)

for(i in 1:length(auc_infos)) {
  p <- p + annotate("text", x = 0.55, y = 0.1 + (i-1)*0.05, label = auc_infos[i], size = 3, hjust = 0, color = "black")
}

print(p)
>>>>>>> 027b93a2ffda7f6e9e9524ea22593a1cf937035a
