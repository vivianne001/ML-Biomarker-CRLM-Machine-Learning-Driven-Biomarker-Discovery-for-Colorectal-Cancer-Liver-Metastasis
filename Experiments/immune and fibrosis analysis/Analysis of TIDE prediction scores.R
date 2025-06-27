################################################################################
# Differential Analysis of TIDE Prediction Scores
################################################################################

# Section 1: Setup and Environment Preparation

rm(list = ls())

# Load required libraries
# install.packages(c("ggplot2", "ggsignif", "dplyr")) # Uncomment if needed
library(ggplot2)
library(ggsignif)
library(dplyr)

# Section 2: Data Loading and Preparation

# Load TIDE prediction results
tide <- read.csv("ACMSD high and low TIDE.csv")

# This script assumes that a phenotype dataframe 'ph' and a 'grouplist' vector
# are already loaded in the environment. 'grouplist' should contain the "high" or "low" classification for each sample.
# Example:
# load("phenotype_data.Rdata")
# ph$group <- grouplist

# Merge TIDE scores with group information
common_samples <- intersect(rownames(ph), tide$Patient)
ph_subset <- ph[common_samples, ]
tide_subset <- tide %>% filter(Patient %in% common_samples)
tide_subset$group <- ph_subset$group

# Set factor levels for plotting
tide_subset$group <- factor(tide_subset$group, levels = c("low", "high"))

# Section 3: Visualization of TIDE Metrics

# Define color palette for groups
group_colors <- c("low" = "#4979b6", "high" = "#d9352a")

# Plot 1: TIDE Score Comparison
plot_tide <- ggplot(tide_subset, aes(x = group, y = TIDE, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1.5) +
  scale_fill_manual(values = group_colors) +
  geom_signif(
    comparisons = list(c("low", "high")),
    test = "wilcox.test",
    map_signif_level = TRUE,
    y_position = max(tide_subset$TIDE, na.rm = TRUE) * 1.05
  ) +
  labs(
    title = "TIDE Score by Expression Group",
    x = "Expression Group",
    y = "TIDE Score"
  ) +
  theme_classic() +
  theme(legend.position = "none")

print(plot_tide)

# Plot 2: Dysfunction Score Comparison
plot_dysfunction <- ggplot(tide_subset, aes(x = group, y = Dysfunction, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1.5) +
  scale_fill_manual(values = group_colors) +
  geom_signif(
    comparisons = list(c("low", "high")),
    test = "wilcox.test",
    map_signif_level = TRUE,
    y_position = max(tide_subset$Dysfunction, na.rm = TRUE) * 1.05
  ) +
  labs(
    title = "Dysfunction Score by Expression Group",
    x = "Expression Group",
    y = "Dysfunction Score"
  ) +
  theme_classic() +
  theme(legend.position = "none")

print(plot_dysfunction)