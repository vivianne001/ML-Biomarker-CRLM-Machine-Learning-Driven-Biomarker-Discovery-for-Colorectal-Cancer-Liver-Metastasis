######################## Comparison of AUC Values for Feature Variables ###################### 
## GSE50760
dat = as.data.frame(t(exp))
dat$group = ph$source_name_ch1
dat$group = ifelse(dat$group == "metastasized cancer", "metastasis", "primary")
table(dat$group)

# Load required packages
library(ROCR)
library(pROC)
library(ggplot2)

# Prepare grouping labels for ROC analysis
labels <- as.numeric(dat$group == "metastasis")

# Compute ROC curves and AUC values for each biomarker
pred_acmsd <- prediction(dat$ACMSD, labels)
perf_acmsd <- performance(pred_acmsd, "tpr", "fpr")
auc_acmsd <- performance(pred_acmsd, "auc")@y.values[[1]]
ci_acmsd <- ci.auc(roc(labels, dat$ACMSD))

pred_serpina3 <- prediction(dat$SERPINA3, labels)
perf_serpina3 <- performance(pred_serpina3, "tpr", "fpr")
auc_serpina3 <- performance(pred_serpina3, "auc")@y.values[[1]]
ci_serpina3 <- ci.auc(roc(labels, dat$SERPINA3))

pred_ptgis <- prediction(dat$PTGIS, labels)
perf_ptgis <- performance(pred_ptgis, "tpr", "fpr")
auc_ptgis <- performance(pred_ptgis, "auc")@y.values[[1]]
ci_ptgis <- ci.auc(roc(labels, dat$PTGIS))

pred_slc2a14 <- prediction(dat$SLC2A14, labels)
perf_slc2a14 <- performance(pred_slc2a14, "tpr", "fpr")
auc_slc2a14 <- performance(pred_slc2a14, "auc")@y.values[[1]]
ci_slc2a14 <- ci.auc(roc(labels, dat$SLC2A14))

# Create separate data frames for each ROC curve
roc_data_acmsd <- data.frame(
  fpr = perf_acmsd@x.values[[1]],
  tpr = perf_acmsd@y.values[[1]],
  group = "ACMSD"
)

roc_data_serpina3 <- data.frame(
  fpr = perf_serpina3@x.values[[1]],
  tpr = perf_serpina3@y.values[[1]],
  group = "SERPINA3"
)

roc_data_ptgis <- data.frame(
  fpr = perf_ptgis@x.values[[1]],
  tpr = perf_ptgis@y.values[[1]],
  group = "PTGIS"
)

roc_data_slc2a14 <- data.frame(
  fpr = perf_slc2a14@x.values[[1]],
  tpr = perf_slc2a14@y.values[[1]],
  group = "SLC2A14"
)

# Combine all ROC curve data
roc_data <- rbind(roc_data_acmsd, roc_data_serpina3, roc_data_ptgis, roc_data_slc2a14)

# Plot ROC curves for all candidate biomarkers
ggplot(roc_data, aes(x = fpr, y = tpr, color = group)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "ROC Curve", x = "1 - Specificity", y = "Sensitivity") +
  scale_color_manual(values = c("ACMSD" = "#EF7F51", "SERPINA3" = "#78D3AC", "PTGIS" = "#9355B0", "SLC2A14" = "#74C1F0")) +
  theme(legend.title = element_blank()) +
  annotate("text", x = 0.6, y = 0.4, label = paste0("ACMSD AUC: ", round(auc_acmsd, 2), " [", round(ci_acmsd[1], 3), "-", round(ci_acmsd[3], 3), "]")) +
  annotate("text", x = 0.6, y = 0.35, label = paste0("SERPINA3 AUC: ", round(auc_serpina3, 2), " [", round(ci_serpina3[1], 3), "-", round(ci_serpina3[3], 3), "]")) +
  annotate("text", x = 0.6, y = 0.3, label = paste0("PTGIS AUC: ", round(auc_ptgis, 2), " [", round(ci_ptgis[1], 3), "-", round(ci_ptgis[3], 3), "]")) +
  annotate("text", x = 0.6, y = 0.25, label = paste0("SLC2A14 AUC: ", round(auc_slc2a14, 2), " [", round(ci_slc2a14[1], 3), "-", round(ci_slc2a14[3], 3), "]"))
