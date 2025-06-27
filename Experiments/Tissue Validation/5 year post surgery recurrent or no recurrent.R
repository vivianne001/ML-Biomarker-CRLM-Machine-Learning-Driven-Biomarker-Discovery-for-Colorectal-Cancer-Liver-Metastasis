###################### Visualization of IHC Score in Recurrent and Non-Recurrent Groups #########################

rm(list = ls())

library(readxl)
library(dplyr)
library(ggsignif)
library(ggplot2)

# Load IHC score data
exp = read_xlsx("Supplementary Tables.xlsx", sheet = 11)
exp = as.data.frame(exp)
exp = exp[, -1]

colnames(exp) = exp[1,]
exp = exp[-1,]
exp = as.data.frame(exp)

# Convert IHC Score to numeric
exp$`IHC Score` = as.numeric(exp$`IHC Score`)

# Prepare data
dat = exp

# Wilcoxon test for IHC Score by Recurrence
wilcox_result <- wilcox.test(dat$`IHC Score`[dat$Recurrence == "Positive"], 
                             dat$`IHC Score`[dat$Recurrence == "Negative"])

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

# Ensure Recurrence is a factor with correct order
dat$Recurrence <- factor(dat$Recurrence, levels = c("Positive", "Negative"))

# Count samples per group
group_counts <- dat %>%
  group_by(Recurrence) %>%
  summarise(count = n())

# Boxplot for IHC Score by Recurrence group
ggplot(dat, aes(x = Recurrence, y = `IHC Score`, fill = Recurrence)) +
  geom_boxplot(width = 0.75, outlier.shape = NA) +
  scale_fill_manual(values = c("Negative" = "red", "Positive" = "lightblue")) +
  geom_signif(comparisons = list(c("Negative", "Positive")), 
              test = "wilcox.test", 
              map_signif_level = FALSE) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 14)
  ) +
  ylab("ACMSD IHC Score") +
  xlab("Group: Recurrent & Non Recurrent")