###################### Heatmap Analysis of Fibrosis and TGFB Pathway Markers  #########################
# GSE204805
rm(list = ls())

# Load marker gene sets and correlation results
load("TGFB and Fibrosis markers.rdata")
load("cor of ACMSD in CRLM samples.rdata")

# Select genes of interest
genes <- c(TGFB, fibrosis)
TGFB <- c("BMP10", "BMP5", "SMAD5", "TGFBR2", "SMAD4", "BMP6", "TGFB1", "SMAD3")
results <- results[results$column %in% TGFB, ]

# Assign significance symbols based on p-values
splitp_vec <- as.vector(results$p.value)
mark_vec <- case_when(
  splitp_vec < 0.001 ~ "***",
  splitp_vec < 0.01 ~ "**",
  splitp_vec < 0.05 ~ "*",
  TRUE ~ ""
)

# Prepare data for plotting
library(ggplot2)
library(reshape2)

results$significance <- mark_vec

df_plot <- data.frame(
  Gene = results$column,
  Correlation = results$correlation,
  Significance = results$significance
)

# Sort genes by correlation (optional)
df_plot <- df_plot[order(df_plot$Correlation, decreasing = FALSE), ]
df_plot$Gene <- factor(df_plot$Gene, levels = df_plot$Gene)

# Draw heatmap
ggplot(df_plot, aes(x = 1, y = Gene, fill = Correlation)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = Significance), size = 5) +
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "firebrick3",
    midpoint = 0.1,
    limits = c(0.1, 0.8),
    name = "Correlation"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Gene Correlation Heatmap with Significance") +
  coord_fixed()