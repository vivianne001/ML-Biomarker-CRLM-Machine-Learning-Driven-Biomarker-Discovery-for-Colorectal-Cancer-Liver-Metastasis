<<<<<<< HEAD
# Visualization of Baseline ACC and AUC Results for Multiple Machine Learning Models

rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(limma)
library(RColorBrewer)
library(readxl)

# Import the result data
result <- read_xlsx("multi learners classification result.xlsx", sheet = 2, col_names = TRUE)
result <- as.data.frame(result)

# Prepare heatmap data
heatmap_data <- result[, c(2, 3, 4)]
rownames(heatmap_data) <- heatmap_data$Name
heatmap_data <- heatmap_data[, -1]
colnames(heatmap_data) <- c("ACC", "AUC")

# Define color palette for the heatmap
col_fun <- colorRamp2(
  c(0.5, 0.75, 1),
  c("#58aaa1", "#FFFFFF", "#e7bc6a")
)
color_vector <- c("#FB8072", "#B3DE69")

# Create column annotation
col_annotation <- columnAnnotation(
  Cohort = factor(colnames(heatmap_data)),
  col = list(
    Cohort = c("AUC" = color_vector[2], "ACC" = color_vector[1])
  ),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 13, fontfamily = "serif", fontface = "bold"),
    labels_gp = gpar(fontsize = 10, fontfamily = "serif"),
    grid_height = unit(0.6, "cm"), grid_width = unit(0.6, "cm")
  ),
  border = FALSE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 13, fontfamily = "serif", fontface = "bold")
)

# Order by AUC for heatmap visualization
heatmap_data <- heatmap_data[order(-heatmap_data$AUC), ]
mean_sort <- sort(heatmap_data$AUC, decreasing = TRUE)

# Create row annotation with barplot
row_bar <- rowAnnotation(
  bar = anno_barplot(mean_sort, bar_width = 0.8, border = FALSE,
                     gp = gpar(fill = "#aa9b81", col = NA),
                     add_numbers = FALSE,
                     width = unit(2, "cm")
  ),
  show_annotation_name = FALSE
)

# Save heatmap as PDF for the validation set
pdf("machine_complexheatmap-1-validationset.pdf", width = 10, height = 20)
Heatmap(as.matrix(heatmap_data),
        col = col_fun,
        show_column_names = FALSE,
        top_annotation = col_annotation,
        right_annotation = row_bar,
        row_title = NULL,
        column_title = NULL,
        column_split = factor(colnames(heatmap_data), levels = colnames(heatmap_data)),
        row_split = factor(rownames(heatmap_data), levels = rownames(heatmap_data)),
        column_gap = unit(2, "mm"),
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 17, fontfamily = "serif", fontface = "bold"),
        row_title_gp = gpar(fontsize = 17, fontfamily = "serif", fontface = "bold"),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 12, fontfamily = "serif", fontface = "bold"),
        heatmap_legend_param = list(
          title = "AUC",
          title_gp = gpar(fontsize = 13, fontfamily = "serif", fontface = "bold"),
          labels_gp = gpar(fontsize = 10, fontfamily = "serif", fontface = "bold"),
          legend_height = unit(8, "cm"),
          grid_width = unit(0.8, "cm")
        ),
        border = FALSE,
        cell_fun = function(j, i, x, y, w, h, col) {
          grid.text(label = format(heatmap_data[i, j], digits = 3, nsmall = 4),
                    x, y, gp = gpar(fontsize = 10))
        }
)
dev.off()

# Basic annotation for columns
annotation_col <- data.frame(
  Cohort = colnames(heatmap_data)
)
ann_colors <- list(
  Cohort = c("AUC" = color_vector[2], "ACC" = color_vector[1])
)

# Save simple heatmap as PDF for the test set
pdf("machine_heatmap of testset.pdf", width = 8, height = 20)
pheatmap(
  heatmap_data,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  display_numbers = TRUE,
  number_format = "%.4f",
  fontsize_row = 11,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  heatmap_legend_param = list(
    title = "AUC",
    title_gp = gpar(fontsize = 13),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(6, "cm"),
    grid_width = unit(0.5, "cm")
  )
)
=======
# Visualization of Baseline ACC and AUC Results for Multiple Machine Learning Models

rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(limma)
library(RColorBrewer)
library(readxl)

# Import the result data
result <- read_xlsx("multi learners classification result.xlsx", sheet = 2, col_names = TRUE)
result <- as.data.frame(result)

# Prepare heatmap data
heatmap_data <- result[, c(2, 3, 4)]
rownames(heatmap_data) <- heatmap_data$Name
heatmap_data <- heatmap_data[, -1]
colnames(heatmap_data) <- c("ACC", "AUC")

# Define color palette for the heatmap
col_fun <- colorRamp2(
  c(0.5, 0.75, 1),
  c("#58aaa1", "#FFFFFF", "#e7bc6a")
)
color_vector <- c("#FB8072", "#B3DE69")

# Create column annotation
col_annotation <- columnAnnotation(
  Cohort = factor(colnames(heatmap_data)),
  col = list(
    Cohort = c("AUC" = color_vector[2], "ACC" = color_vector[1])
  ),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 13, fontfamily = "serif", fontface = "bold"),
    labels_gp = gpar(fontsize = 10, fontfamily = "serif"),
    grid_height = unit(0.6, "cm"), grid_width = unit(0.6, "cm")
  ),
  border = FALSE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 13, fontfamily = "serif", fontface = "bold")
)

# Order by AUC for heatmap visualization
heatmap_data <- heatmap_data[order(-heatmap_data$AUC), ]
mean_sort <- sort(heatmap_data$AUC, decreasing = TRUE)

# Create row annotation with barplot
row_bar <- rowAnnotation(
  bar = anno_barplot(mean_sort, bar_width = 0.8, border = FALSE,
                     gp = gpar(fill = "#aa9b81", col = NA),
                     add_numbers = FALSE,
                     width = unit(2, "cm")
  ),
  show_annotation_name = FALSE
)

# Save heatmap as PDF for the validation set
pdf("machine_complexheatmap-1-validationset.pdf", width = 10, height = 20)
Heatmap(as.matrix(heatmap_data),
        col = col_fun,
        show_column_names = FALSE,
        top_annotation = col_annotation,
        right_annotation = row_bar,
        row_title = NULL,
        column_title = NULL,
        column_split = factor(colnames(heatmap_data), levels = colnames(heatmap_data)),
        row_split = factor(rownames(heatmap_data), levels = rownames(heatmap_data)),
        column_gap = unit(2, "mm"),
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 17, fontfamily = "serif", fontface = "bold"),
        row_title_gp = gpar(fontsize = 17, fontfamily = "serif", fontface = "bold"),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 12, fontfamily = "serif", fontface = "bold"),
        heatmap_legend_param = list(
          title = "AUC",
          title_gp = gpar(fontsize = 13, fontfamily = "serif", fontface = "bold"),
          labels_gp = gpar(fontsize = 10, fontfamily = "serif", fontface = "bold"),
          legend_height = unit(8, "cm"),
          grid_width = unit(0.8, "cm")
        ),
        border = FALSE,
        cell_fun = function(j, i, x, y, w, h, col) {
          grid.text(label = format(heatmap_data[i, j], digits = 3, nsmall = 4),
                    x, y, gp = gpar(fontsize = 10))
        }
)
dev.off()

# Basic annotation for columns
annotation_col <- data.frame(
  Cohort = colnames(heatmap_data)
)
ann_colors <- list(
  Cohort = c("AUC" = color_vector[2], "ACC" = color_vector[1])
)

# Save simple heatmap as PDF for the test set
pdf("machine_heatmap of testset.pdf", width = 8, height = 20)
pheatmap(
  heatmap_data,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  display_numbers = TRUE,
  number_format = "%.4f",
  fontsize_row = 11,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  heatmap_legend_param = list(
    title = "AUC",
    title_gp = gpar(fontsize = 13),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(6, "cm"),
    grid_width = unit(0.5, "cm")
  )
)
>>>>>>> 027b93a2ffda7f6e9e9524ea22593a1cf937035a
dev.off()