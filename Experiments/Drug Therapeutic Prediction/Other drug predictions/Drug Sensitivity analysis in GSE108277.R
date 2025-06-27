###################### Visualization of Drug Sensitivity Differences in GSE108277 #########################
# GSE108277
rm(list = ls())

library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggsignif)
library(ggplot2)

# Log2 transformation of expression data
# load("GSE108277_processed.Rdata")
exp = log2(exp + 1)
range(exp)

# Filter for treatment and response groups
table(ph$`treatment:ch1`)
ph = ph[ph$`treatment:ch1` %in% c("acute", "chronic"), ]

table(ph$`cetuximab response:ch1`)
ph = ph[ph$`cetuximab response:ch1` %in% c("PR", "SD"), ]

# Match samples between expression and phenotype data
k = intersect(rownames(ph), colnames(exp))
exp = exp[, k]
ph = ph[k, ]
identical(rownames(ph), colnames(exp))

# Assign group labels
ph$group = ph$`cetuximab response:ch1`
ph$group = ifelse(ph$group == "PR", "Partial_Response", "Stable_Disease")

# Prepare data for boxplot
dat = t(exp)
dat = as.data.frame(dat)
dat$group = ph$group
dat$group = as.factor(dat$group)

# t-test and Wilcoxon test for ACMSD expression
t_test_result <- t.test(dat$ACMSD[dat$group == "Partial_Response"], 
                        dat$ACMSD[dat$group == "Stable_Disease"])
wilcox_result <- wilcox.test(dat$ACMSD[dat$group == "Partial_Response"], 
                             dat$ACMSD[dat$group == "Stable_Disease"])

# Function to map p-value to significance symbol
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
significance_symbols <- map_pvalue_to_signif_mark(t_test_result$p.value)
significance_symbols <- map_pvalue_to_signif_mark(wilcox_result$p.value)

# Ensure group factor order
dat$group <- factor(dat$group, levels = c("Partial_Response", "Stable_Disease"))

# Boxplot for ACMSD expression by response group
ggplot(dat, aes(x = group, y = ACMSD, fill = group)) +
  geom_boxplot(width = 0.75, outlier.shape = NA) +
  scale_fill_manual(values = c("Stable_Disease" = "red", "Partial_Response" = "lightblue")) +
  geom_signif(comparisons = list(c("Stable_Disease", "Partial_Response")), 
              test = "t.test", 
              map_signif_level = TRUE) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14)) +
  ylab("ACMSD expression (log2)")

# Prepare data for multi-gene boxplot
dat = t(exp)
dat = as.data.frame(dat)
dat$group = ph$`cetuximab response:ch1`
dat = dat[, c("SERPINA3", "PTGIS", "SLC2A14", "ACMSD", "group")]
