rm(list = ls())

###################### Selection of candidate biomarkers ########################
# GSE204805
# Transform prognosis information to numeric and logical format
table(ph$characteristics_ch1.2)
ph = ph[, c("geo_accession", "characteristics_ch1.2")]
ph$group = ifelse(ph$characteristics_ch1.2 == "cell type: Metastasis", "metastasis", "Primary")
table(ph$group)

# Log2 transformation of expression data and preprocessing
exp = log2(exp + 1)
exp = t(exp)
exp = as.data.frame(exp)
exp = exp[, c("SERPINA3", "PTGIS", "ACMSD", "SLC2A14")]
dat = exp

# Select common samples for pairing phenotype and expression data
k = intersect(rownames(dat), rownames(ph))
ph = ph[k, ]
dat = dat[k, ]
dat$group = ph$group

head(dat)

#############################################
# Ensure PrimaryOC package is installed and load pROC
if (!require(PrimaryOC)) {
  install.packages("PrimaryOC")
}
library(pROC)

# Prepare data for ROC analysis
mpe = t(exp)
mpe = as.data.frame(mpe)
k = intersect(rownames(mpe), rownames(ph))
mpe = mpe[k, ]
ph = ph[k, ]
table(ph$`cell type:ch1`)
mpe$group = ph$`cell type:ch1`
mpe$group = ifelse(mpe$group == "LM", "Metastasis", "Primary")
mpe$group <- factor(mpe$group, levels = c("Metastasis", "Primary"))

# Calculate and plot ROC curves for selected genes
for (gene in c("ACMSD", "GPA33", "CDH17", "SATB2", "KRT7", "KRT15", "KRT20", "KRT18", "TERT")) {
  roc_obj <- roc(mpe$group, mpe[[gene]], levels=c("Metastasis","Primary"))
  plot(roc_obj, main=paste("ROC curve for", gene))
  cat("AUC for", gene, ":", auc(roc_obj), "\n")
}
for (gene in c("ACMSD", "KRAS", "NRAS", "BRAF", "PIK3CA", "SMAD4", "APC", "PTEN")) {
  roc_obj <- roc(mpe$group, mpe[[gene]], levels=c("Metastasis","Primary"))
  plot(roc_obj, main=paste("ROC curve for", gene))
  cat("AUC for", gene, ":", auc(roc_obj), "\n")
}
for (gene in c("ACMSD", "SERPINA3", "PTGIS", "SLC2A14")) {
  roc_obj <- roc(mpe$group, mpe[[gene]], levels=c("Metastasis","Primary"))
  plot(roc_obj, main=paste("ROC curve for", gene))
  cat("AUC for", gene, ":", auc(roc_obj), "\n")
}

###########################
library(PrimaryOC)
# Create empty data frame to store AUC results
auc_results <- data.frame(Gene = character(), AUC = numeric(), stringsAsFactors = FALSE)
# Calculate ROC and AUC for each gene and store the results
for (gene in c("SERPINA3", "PTGIS", "ACMSD", "SLC2A14")) {
  roc_obj <- roc(mpe$group, mpe[[gene]], levels=c("Metastasis", "Primary"))
  auc_value <- auc(roc_obj)
  auc_results <- rbind(auc_results, data.frame(Gene = gene, AUC = auc_value))
}
# Print the results
print(auc_results)

#######################################################
# Additional ROC curve analysis for a single biomarker (ACMSD)
library(pROC)
mpe = t(exp)
mpe = as.data.frame(mpe)
k = intersect(rownames(mpe), rownames(ph))
mpe = mpe[k, ]
ph = ph[k, ]
mpe$group = ph$`metastatic tumor site:ch1`
mpe$group <- factor(mpe$group, levels = c("Metastasis", "Primary"))
roc_obj <- roc(mpe$group, mpe[["ACMSD"]], levels=c("Metastasis", "Primary"))
plot(roc_obj, main="ROC curve for ACMSD")
cat("AUC for ACMSD:", auc(roc_obj), "\n")

#####################################################
## Single biomarker ROC curve visualization
library(pROC)
library(ggplot2)

# Calculate ROC data for ACMSD
roc_data <- roc(mpe$group, mpe[["ACMSD"]], levels=c("Metastasis","Primary"))
roc_df <- data.frame(
  specificity = roc_data$specificities,
  sensitivity = roc_data$sensitivities
)
auc_value <- auc(roc_data)
auc_ci <- ci.auc(roc_data)
auc_text <- sprintf("AUC = %.2f (95%% CI: %.2f-%.2f)", auc_value, auc_ci[1], auc_ci[3])
print(auc_text)

# Plot ROC using ggplot2
p <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line() +
  geom_abline(linetype = "dashed", color = "grey") +
  scale_x_continuous(limits = c(0, 1), name = "1-Specificity") +
  scale_y_continuous(limits = c(0, 1), name = "Sensitivity") +
  ggtitle("ROC curve for ACMSD") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black")) +
  annotate("text", x = 0.6, y = 0.15, label = "mCRC Liver Metastasis & Primary", size = 4, hjust = 0) +
  annotate("text", x = 0.6, y = 0.1, label = auc_text, size = 4, hjust = 0)
print(p)

####################################
# ROC curves for multiple genes in a combined plot
library(pROC)
library(ggplot2)
genes <- c("ACMSD", "SERPINA3", "PTGIS", "SLC2A14")
colors <- c("SERPINA3" = "#ADAFB1", "PTGIS" = "#E08D8B", "SLC2A14" = "#ECB884", "ACMSD" = "#9584C1")
roc_df <- data.frame()
auc_infos <- c()
for(gene in genes) {
  roc_data <- roc(mpe$group, mpe[[gene]], levels=c("Metastasis", "Primary"))
  temp_df <- data.frame(
    specificity = roc_data$specificities,
    sensitivity = roc_data$sensitivities,
    gene = gene
  )
  roc_df <- rbind(roc_df, temp_df)
  auc_value <- auc(roc_data)
  ci <- ci.auc(roc_data)
  auc_info <- sprintf("%s: AUC = %.2f (95%% CI: %.2f-%.2f)", gene, auc_value, ci[1], ci[3])
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
  p <- p + annotate("text", x = 0.65, y = 0.1 + (i-1)*0.05, label = auc_infos[i], size = 3, hjust = 0, color = "black")
}
print(p)