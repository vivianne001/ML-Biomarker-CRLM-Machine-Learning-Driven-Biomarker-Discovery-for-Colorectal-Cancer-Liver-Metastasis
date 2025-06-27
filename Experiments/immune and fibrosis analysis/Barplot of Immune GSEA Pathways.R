######################## Barplot of Immune GSEA Pathways ########################

library(readxl)
library(dplyr)
library(ggplot2)

# Load pathway enrichment results
pathway <- read_xlsx("Immune_Inflammatory_GSEA_Statistics.xlsx", sheet = 1, col_names = TRUE)
pathway <- as.data.frame(pathway)

# Add significance symbols and color coding
pathway <- pathway %>%
  mutate(
    significance = case_when(
      pvalue < 0.001 ~ "***",
      pvalue < 0.01 ~ "**",
      pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    color = ifelse(base == "GO", "lightcoral", "darkred")
  )

# Create barplot
ggplot(pathway, aes(x = factor(Pathway, levels = Pathway), y = NES, fill = color)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = significance), vjust = -0.5, size = 5) +
  scale_fill_identity() +
  labs(
    x = "Pathway",
    y = "NES",
    title = "Barplot of Pathways with Significance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  )