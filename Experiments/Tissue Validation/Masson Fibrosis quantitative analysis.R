######################### Quantitative Analysis and Visualization of Masson Fibrosis Score ############################

rm(list = ls())

library(readxl)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(dplyr)

# Load data
exp = read_xlsx("Supplementary Tables.xlsx", sheet = 16)
exp = as.data.frame(exp)
colnames(exp) = exp[1, ]
exp = exp[-1, ]
exp = exp[, -1]

# Convert Masson Area and IHC Score to numeric
exp$`Masson Area` = as.numeric(exp$`Masson Area`)
exp$`Masson Area` = round(exp$`Masson Area`, 3)
exp$`IHC Score` = as.numeric(exp$`IHC Score`)

# Assign group2 based on IHC Score median
exp$group2 = ifelse(exp$`IHC Score` > median(exp$`IHC Score`), "high", "low")
exp$group2 = as.factor(exp$group2)

# Convert IHC Score to categorical labels for group
exp$`IHC Score` = gsub("0", "Negative", exp$`IHC Score`)
exp$`IHC Score` = gsub("1", "Weak Positive", exp$`IHC Score`)
exp$`IHC Score` = gsub("2", "Positive", exp$`IHC Score`)
exp$`IHC Score` = gsub("4|6", "Strong Positive", exp$`IHC Score`)

# Prepare data for plotting
dat = exp
dat$group = exp$`IHC Score`
colnames(dat)[12] = "masson"
dat$group <- factor(dat$group, levels = c("Negative", "Weak Positive", "Positive", "Strong Positive"))

# Boxplot for Masson Area by IHC Score group
ggplot(dat, aes(group, masson, fill = group)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("Negative", "Weak Positive"),
                                 c("Negative", "Positive"),
                                 c("Negative", "Strong Positive"),
                                 c("Weak Positive", "Positive"),
                                 c("Weak Positive", "Strong Positive"),
                                 c("Positive", "Strong Positive")),
              map_signif_level = TRUE,
              textsize = 4.5,
              test = wilcox.test,
              step_increase = 0.1) +
  scale_fill_manual(values = c("Negative" = "#78D3AC", "Weak Positive" = "#EF7F51",
                               "Positive" = "#9355B0", "Strong Positive" = "darkred")) +
  guides(fill = FALSE) +
  xlab(NULL) +
  ylab("Fibrosis(%)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))

# Boxplot for Masson Area by IHC Score median group (high/low)
dat$group2 <- factor(dat$group2, levels = c("low", "high"))

ggplot(dat, aes(x = group2, y = masson, fill = group2)) +
  geom_boxplot(width = 0.75, outlier.shape = NA) +
  scale_fill_manual(values = c("high" = "red", "low" = "lightblue")) +
  geom_signif(comparisons = list(c("low", "high")), 
              test = "wilcox.test", 
              map_signif_level = TRUE) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14)) +
  ylab("Fibrosis(%)") +
  xlab("Group: ACMSD IHC Protein Level")

# Statistical tests
# t_test_result <- t.test(dat$masson[dat$group2 == "high"], dat$masson[dat$group2 == "low"])
wilcox_result <- wilcox.test(dat$masson[dat$group2 == "high"], dat$masson[dat$group2 == "low"])

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

# significance_symbols <- map_pvalue_to_signif_mark(t_test_result$p.value)
significance_symbols <- map_pvalue_to_signif_mark(wilcox_result$p.value)