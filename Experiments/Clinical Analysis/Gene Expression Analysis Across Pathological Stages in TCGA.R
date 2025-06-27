
############ Gene Expression Analysis Across Pathological Stages in TCGA ##############
################################################################################

# Section 1: Setup and Environment Preparation

rm(list = ls())

# Load required libraries
# install.packages(c("tidyverse", "ggsignif", "patchwork")) # Uncomment if needed
library(tidyverse)
library(ggsignif)
library(patchwork)

# Section 2: Data Loading and Initial Cleaning

# This script assumes 'exp' (expression matrix) and 'clin' (clinical data) are loaded.
# Example: load("TCGA_COADREAD_Raw Counts.Rdata")

exp <- as.data.frame(exp)
clin <- as.data.frame(clin)

# Log2 transform expression data
exp <- log2(exp + 1)

# Ensure all expression data is numeric
exp <- data.frame(lapply(exp, function(x) as.numeric(as.character(x))),
                  row.names = rownames(exp))
exp[is.na(exp)] <- 0

# Standardize sample names
colnames(exp) <- gsub("\\.", "-", colnames(exp))
colnames(exp) <- sub("-\\d{2}$", "", colnames(exp))
rownames(clin) <- clin$PATIENT_ID

# Section 3: Data Harmonization and Merging

common_samples <- intersect(colnames(exp), rownames(clin))
exp_subset <- exp[, common_samples]
clin_subset <- clin[common_samples, ]

# Ensure sample order matches
identical(colnames(exp_subset), rownames(clin_subset))

# Transpose expression data: samples as rows, genes as columns
dat <- as.data.frame(t(exp_subset))

# Section 4: Clinical Variable Processing

# AJCC Pathologic Stage grouping
clin_stage <- clin_subset
clin_stage$AJCC_PATHOLOGIC_TUMOR_STAGE <- gsub("STAGE IA|STAGE IB", "STAGE I", clin_stage$AJCC_PATHOLOGIC_TUMOR_STAGE)
clin_stage$AJCC_PATHOLOGIC_TUMOR_STAGE <- gsub("STAGE IIA|STAGE IIB|STAGE IIC", "STAGE II", clin_stage$AJCC_PATHOLOGIC_TUMOR_STAGE)
clin_stage$AJCC_PATHOLOGIC_TUMOR_STAGE <- gsub("STAGE IIIA|STAGE IIIB|STAGE IIIC", "STAGE III", clin_stage$AJCC_PATHOLOGIC_TUMOR_STAGE)
clin_stage$AJCC_PATHOLOGIC_TUMOR_STAGE <- gsub("STAGE IVA|STAGE IVB", "STAGE IV", clin_stage$AJCC_PATHOLOGIC_TUMOR_STAGE)
clin_stage$AJCC_PATHOLOGIC_TUMOR_STAGE <- gsub("STAGE IV|STAGE III", "STAGE III/IV", clin_stage$AJCC_PATHOLOGIC_TUMOR_STAGE)
clin_stage <- clin_stage[!grepl("\\[|\\]", clin_stage$AJCC_PATHOLOGIC_TUMOR_STAGE), ]

# TNM Pathologic T-Stage grouping
clin_t_stage <- clin_subset
clin_t_stage$PATH_T_STAGE <- gsub("T4A|T4B", "T4", clin_t_stage$PATH_T_STAGE)
clin_t_stage$PATH_T_STAGE <- gsub("T1|T2", "T1/2", clin_t_stage$PATH_T_STAGE)
clin_t_stage <- clin_t_stage[clin_t_stage$PATH_T_STAGE != "TIS", ]

# TNM Pathologic N-Stage grouping
clin_n_stage <- clin_subset
clin_n_stage$PATH_N_STAGE <- gsub("N1A|N1B|N1C", "N1", clin_n_stage$PATH_N_STAGE)
clin_n_stage$PATH_N_STAGE <- gsub("N2A|N2B|NX", "N2/X", clin_n_stage$PATH_N_STAGE)

# Section 5: Visualization

# Reusable plotting function
create_expression_plot <- function(data, gene_name, group_name, comparisons, colors, y_label, plot_type = "boxplot") {
  p <- ggplot(data, aes_string(x = group_name, y = gene_name, fill = group_name))
  if (plot_type == "violin") {
    p <- p + geom_violin()
  } else {
    p <- p + geom_boxplot()
  }
  p <- p +
    geom_signif(
      comparisons = comparisons,
      map_signif_level = TRUE,
      test = "t.test",
      step_increase = 0.1,
      textsize = 4.5
    ) +
    scale_fill_manual(values = colors) +
    guides(fill = "none") +
    labs(x = NULL, y = y_label) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
  return(p)
}

# Plot 1: Gene Expression by AJCC Stage (Boxplot)
dat_stage <- dat[rownames(clin_stage), ]
dat_stage$group <- clin_stage$AJCC_PATHOLOGIC_TUMOR_STAGE
dat_stage <- dat_stage[, c("SERPINA3", "PTGIS", "ACMSD", "SLC2A14", "group")]
dat_stage <- na.omit(dat_stage)

plot1 <- create_expression_plot(
  data = dat_stage,
  gene_name = "ACMSD",
  group_name = "group",
  comparisons = list(c("STAGE I", "STAGE II"), c("STAGE I", "STAGE III/IV")),
  colors = c("STAGE I" = "#EF7F51", "STAGE II" = "#78D3AC", "STAGE III/IV" = "#9355B0"),
  y_label = "ACMSD Expression (Log2)",
  plot_type = "boxplot"
)
print(plot1)

# Plot 2: Gene Expression by T-Stage (Violin Plot)
dat_t_stage <- dat[rownames(clin_t_stage), ]
dat_t_stage$group <- clin_t_stage$PATH_T_STAGE
dat_t_stage <- dat_t_stage[, c("SERPINA3", "PTGIS", "ACMSD", "SLC2A14", "group")]
dat_t_stage <- na.omit(dat_t_stage)

plot2 <- create_expression_plot(
  data = dat_t_stage,
  gene_name = "ACMSD",
  group_name = "group",
  comparisons = list(c("T1/2", "T3"), c("T1/2", "T4")),
  colors = c("T1/2" = "#EF7F51", "T3" = "#78D3AC", "T4" = "#9355B0"),
  y_label = "ACMSD Expression (Log2)",
  plot_type = "violin"
)
print(plot2)