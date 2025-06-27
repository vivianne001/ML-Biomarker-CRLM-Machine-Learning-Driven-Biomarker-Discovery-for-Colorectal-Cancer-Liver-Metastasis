
########## Correlation Analysis of Gene Expression with Tumor Purity ###########
################################################################################

# Section 1: Setup and Environment Preparation

rm(list = ls())

# Load required libraries
# install.packages(c("dplyr", "ggplot2", "ggpubr", "patchwork")) # Uncomment if needed
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

# Section 2: Data Loading and Preparation

# This script assumes that 'tpm' (expression matrix), 'ph' (phenotype data),
# and 'estimate.normalized' (ESTIMATE scores) are loaded.
# load("GSE204805_processed.Rdata")
# load("estimate.tpm.Rdata")
# load("exp of tpm.Rdata")

# Filter for the desired cell type in the phenotype data
ph_filtered <- ph[ph$characteristics_ch1.2 %in% c("cell type: LM"), ]

# Find common samples across all three dataframes
common_samples <- intersect(rownames(ph_filtered), colnames(tpm))
common_samples <- intersect(common_samples, estimate.normalized$ID)

# Subset all data to the common samples
ph_final <- ph_filtered[common_samples, ]
tpm_final <- tpm[, common_samples]
estimate_final <- estimate.normalized %>% filter(ID %in% common_samples)

# Create the final analysis dataframe
dat <- as.data.frame(t(tpm_final)) %>%
  select(SERPINA3, PTGIS, ACMSD, SLC2A14)
dat$TumorPurity <- estimate_final$TumorPurity_estimate

# Section 3: Visualization (Reusable Function)

# Create a reusable plotting function for correlation analysis
create_correlation_plot <- function(data, x_var, y_var, x_lab, y_lab) {
  cor_test <- cor.test(data[[x_var]], data[[y_var]], method = "spearman")
  rho <- round(cor_test$estimate, 2)
  pval <- format.pval(cor_test$p.value, digits = 2, eps = 0.001)
  p <- ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(color = "#0073C2", alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "#E41A1C", fill = "grey80") +
    labs(x = x_lab, y = y_lab) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold")
    ) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = paste0("Spearman's R = ", rho, "\np = ", pval),
      size = 4
    )
  return(p)
}

# Generate plots for each gene
p1 <- create_correlation_plot(dat, "SERPINA3", "TumorPurity", "SERPINA3 Expression", "Tumor Purity Score")
p2 <- create_correlation_plot(dat, "PTGIS", "TumorPurity", "PTGIS Expression", "Tumor Purity Score")
p3 <- create_correlation_plot(dat, "ACMSD", "TumorPurity", "ACMSD Expression", "Tumor Purity Score")
p4 <- create_correlation_plot(dat, "SLC2A14", "TumorPurity", "SLC2A14 Expression", "Tumor Purity Score")

# Combine plots into a single figure
combined_plot <- (p1 | p2) / (p3 | p4)
print(combined_plot)

# Section 4: (Alternative) Automated Loop for Plotting

feature_genes <- c("SERPINA3", "PTGIS", "ACMSD", "SLC2A14")
plot_list <- list()

for (gene in feature_genes) {
  p <- ggplot(dat, aes_string(x = gene, y = "TumorPurity")) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
    labs(x = paste(gene, "Expression"), y = "Tumor Purity Score") +
    theme_minimal() +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top")
  plot_list[[gene]] <- p
}

wrap_plots(plot_list, ncol = 2)