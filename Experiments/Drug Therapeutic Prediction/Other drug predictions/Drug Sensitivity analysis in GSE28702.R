###################### Visualization of Drug Response Proportions and Pie Charts in GSE28702 #########################

rm(list = ls())

library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(scales)

# Log2 transformation of expression data
# load("GSE28702_processed.Rdata")
exp = log2(exp + 1)
range(exp)

# Filter for liver location samples
table(ph$characteristics_ch1.1)
ph = ph[ph$characteristics_ch1.1 %in% c("location: Liver"), ]

# Match samples between expression and phenotype data
k = intersect(rownames(ph), colnames(exp))
ph = ph[k, ]
exp = exp[, k]

# Assign group labels based on treatment
table(ph$`mfolfox6:ch1`)
ph$group = ph$`mfolfox6:ch1`

# Prepare data for plotting
dat = t(exp)
dat = as.data.frame(dat)
dat$group = ph$`cell type:ch1`
dat = dat[, c("SERPINA3", "PTGIS", "SLC2A14", "ACMSD", "group")]

# Dichotomize the first four genes by median
dat[, 1:4] <- apply(dat[, 1:4], 2, function(x) ifelse(x > median(x), "high", "low"))

# Calculate proportions of high/low ACMSD in each group
dat_summary <- dat %>%
  group_by(group, ACMSD) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(prop = count / sum(count))

# Barplot of ACMSD high/low proportions by group
ggplot(dat_summary, aes(x = group, y = prop, fill = ACMSD)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = percent(prop)), position = position_stack(vjust = 0.5)) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    x = "Group",
    y = "Proportion",
    fill = "ACMSD Level",
    title = "Chemotherapy + BEVACIZUMAB Response Proportion Plot of ACMSD"
  ) +
  theme_minimal() +
  theme(axis.line = element_line())

# Pie chart of ACMSD high/low proportions by group
ggplot(dat_summary, aes(x = group, y = prop, fill = ACMSD)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = percent(prop)), position = position_stack(vjust = 0.5)) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    x = NULL,
    y = "Proportion",
    fill = "ACMSD Level",
    title = "FOLFOX6 Response Proportion Plot of ACMSD"
  ) +
  facet_wrap(~ group) +
  theme_minimal() +
  theme(
    axis.line = element_line(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )