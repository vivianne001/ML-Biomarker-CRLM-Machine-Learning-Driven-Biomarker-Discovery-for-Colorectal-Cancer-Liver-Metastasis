rm(list = ls())

###################### Analysis of Intergroup Boxplots #########################
# GSE50760

library("pheatmap")
library("RColorBrewer")
library(dplyr)

# Assign group information to the phenotype dataframe
ph$group = grouplist

# Import candidate gene list
flt = flt[flt$score == 1,]

interest_gene = flt$feature
exp = exp[interest_gene,]
dat = t(exp)
dat = as.data.frame(dat)

# Log2 transformation
dat = log2(dat + 1)
range(dat)
dat$group = ph$group
dat$group = as.factor(dat$group)

# t-test for target gene expression between groups
t_test_result <- t.test(dat$ACMSD[dat$group == "metastasis"],
                        dat$ACMSD[dat$group == "primary"])

# Wilcoxon test for target gene expression between groups
wilcox_result <- wilcox.test(dat$ACMSD[dat$group == "metastasis"],
                             dat$ACMSD[dat$group == "primary"])

# Load necessary libraries for visualization and statistics
library(ggsignif)
library(ggplot2)

# Function to map p-value to significance stars
map_pvalue_to_signif_mark <- function(pvalue) {
  if (pvalue < 0.001) {
    return("***")
  } else if (pvalue < 0.01) {
    return("**")
  } else if (pvalue < 0.055) {
    return("*")
  } else {
    return("ns")
  }
}

# Get significance symbols
significance_symbol_ttest <- map_pvalue_to_signif_mark(t_test_result$p.value)
# significance_symbol_wilcoxon <- map_pvalue_to_signif_mark(wilcox_result$p.value)

# Example: Calculate outliers for a variable (e.g., testPtype.pr$BRD_K30748066)
# outliers <- boxplot.stats(testPtype.pr$BRD_K30748066)$out

# Boxplot for drug sensitivity difference (Single marker)
dat$group <- factor(dat$group, levels = c("primary", "metastasis"))

# Example Boxplot for a specific gene (e.g., SLC2A14)
ggplot(dat, aes(x = group, y = ACMSD, fill = group)) +
  geom_boxplot(width = 0.75, outlier.shape = NA) +
  scale_fill_manual(values = c("metastasis" = "red", "primary" = "lightblue")) +
  geom_signif(comparisons = list(c("metastasis", "primary")),
              test = "t.test",
              map_signif_level = TRUE) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14)) +
  ylab("ACMSD expression (log2+1)") +
  xlab("Group: ACMSD metastasis and primary")


############## Gene Expression Analysis for Multiple Biomarkers ###############
options(stringsAsFactors = FALSE)
library(tidyverse)
library(ggsignif)
library(reshape2)
library(IOBR)

# Reshape data for multi-marker analysis
TME_NEW = melt(dat)
colnames(TME_NEW) = c("Group", "Gene", "Express")
TME_NEW$Express <- as.numeric(TME_NEW$Express)

# Order genes (features) by median expression in the metastasis group
plot_order = TME_NEW[TME_NEW$Group == "metastasis", ] %>%
  group_by(Gene) %>%
  summarise(m = median(Express)) %>%
  arrange(desc(m)) %>%
  pull(Gene)

# Theme for plotting
mytheme <- theme(
  plot.title = element_text(size = 12, color = "black", hjust = 0.5),
  axis.title = element_text(size = 12, color = "black"),
  axis.text = element_text(size = 12, color = "black"),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1),
  panel.grid = element_blank(),
  legend.position = "top",
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 12)
)

# Boxplot for the 29 candidate biomarkers
box_TME <- ggplot(TME_NEW, aes(x = Gene, y = Express)) +
  labs(y = "Expression", x = NULL, title = "JAK-STAT Pathway Mediated Drugs Responsiveness") +
  geom_boxplot(aes(fill = Group), position = position_dodge(0.5), width = 0.5, outlier.alpha = 0) +
  scale_fill_manual(values = c("#1CB4B8", "#EB7369")) +
  theme_classic() +
  mytheme +
  stat_compare_means(aes(group = Group),
                     label = "p.signif",
                     method = "t.test",
                     hide.ns = TRUE) +
  scale_y_continuous(limits = c(-10, 15))

box_TME

# Additional visualization with better panel formatting
library(ggpubr)
TME_NEW$Group <- factor(TME_NEW$Group, levels = c("primary", "metastasis"))

p <- ggplot(TME_NEW, aes(x = Gene, y = Express, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("primary" = "lightblue", "metastasis" = "red"),
                    labels = c("primary" = "primary CRC", "metastasis" = "CRLM")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.line = element_line(color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  stat_compare_means(aes(group = Group), method = "t.test", label = "p.signif", label.y = max(TME_NEW$Express))

print(p)