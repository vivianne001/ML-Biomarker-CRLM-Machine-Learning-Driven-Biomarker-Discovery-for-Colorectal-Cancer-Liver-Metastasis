rm(list = ls())

###################### Intergroup Boxplot Analysis #########################
# GSE204805
library(pheatmap)
library(RColorBrewer)
library(dplyr)

# Assign group information
ph$group = grouplist
table(ph$`cell type:ch1`)
ph = ph[ph$`cell type:ch1` == "LM", ]
k = intersect(rownames(ph), colnames(exp))
exp = exp[, k]
ph = ph[k, ]
identical(rownames(ph), colnames(exp))

# Select candidate genes
flt = flt[flt$score == 1, ]
interest_gene = flt$feature
exp = exp[interest_gene, ]
dat = exp
dat = as.data.frame(dat)

################### Boxplot Analysis ######################

# Prepare data for boxplot
dat = t(dat)
dat = as.data.frame(dat)
dat = log2(dat + 1)
range(dat)
dat$group = ann_col$group
dat$group = as.factor(dat$group)

# t-test and Wilcoxon test of candidate biomarkers for group differences
t_test_result <- t.test(dat$ACMSD[dat$group == "metastasis"],
                        dat$ACMSD[dat$group == "primary"])
wilcox_result <- wilcox.test(dat$ACMSD[dat$group == "metastasis"],
                             dat$ACMSD[dat$group == "primary"])

library(ggsignif)
library(ggplot2)

# Function to map p-value to significance mark
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

# Get significance marks
significance_symbol_ttest <- map_pvalue_to_signif_mark(t_test_result$p.value)
# significance_symbol_wilcoxon <- map_pvalue_to_signif_mark(wilcox_result$p.value)

# Example: Calculate outliers for a variable
# outliers <- boxplot.stats(testPtype.pr$BRD_K30748066)$out

# Boxplot for a single gene (e.g., NADK)
dat$group <- factor(dat$group, levels = c("primary", "metastasis"))
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

################### Multi-Gene Drug Sensitivity Analysis ###################

options(stringsAsFactors = FALSE)
library(tidyverse)
library(reshape2)
library(IOBR)

# Reshape data for multi-gene analysis
TME_NEW = melt(dat)
colnames(TME_NEW) = c("Group", "Gene", "Express")
TME_NEW$Express <- as.numeric(TME_NEW$Express)
head(TME_NEW)

# Order genes by median expression in metastasis group
plot_order = TME_NEW[TME_NEW$Group == "metastasis", ] %>%
  group_by(Gene) %>%
  summarise(m = median(Express)) %>%
  arrange(desc(m)) %>%
  pull(Gene)

# Custom theme for plots
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

# Boxplot for all candidate biomarkers
box_TME <- ggplot(TME_NEW, aes(x = Gene, y = Express)) +
  labs(y = "Expression", x = NULL, title = "JAK-STAT Pathway Mediated Drug Responsiveness") +
  geom_boxplot(aes(fill = Group), position = position_dodge(0.5), width = 0.5, outlier.alpha = 0) +
  scale_fill_manual(values = c("#1CB4B8", "#EB7369")) +
  theme_classic() + mytheme +
  stat_compare_means(aes(group = Group),
                     label = "p.signif",
                     method = "t.test",
                     hide.ns = TRUE) +
  scale_y_continuous(limits = c(-10, 15))
box_TME

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