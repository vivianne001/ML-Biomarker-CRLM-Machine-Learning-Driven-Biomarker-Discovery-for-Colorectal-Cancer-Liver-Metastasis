###################### Identification and Visualization of Functional DEGs #########################

rm(list = ls())

#############################################################################
## Functional DEGs Identification

# Load DEG matrix and metabolic gene list
load("grouplist-volcano-metastasis-DATA.P.Value 0.05.Rdata")
metabolic_gene = read.csv("metabolic genes.csv")

# Extract upregulated and downregulated DEGs
DEGS = DEG[DEG$change %in% c("up", "down"), ]

# Save metabolic gene list as CSV
write.csv(metabolic_gene, file = "metabolic genes.csv")

# Intersect DEGs and metabolic genes to obtain functional DEGs
dd = intersect(DEGS$gene, metabolic_gene$geneSymbols)
function_genes = metabolic_gene[metabolic_gene$geneSymbols %in% dd, ]

# Save functional DEGs as CSV and Rdata
write.csv(function_genes, file = "function dataset.csv")
save(function_genes, file = "functional dataset.Rdata")

##############################################################################
## Venn Plot of Functional DEGs

# Prepare gene names for Venn plot
DEGS_NAME = DEGS$gene_name
meta_genes = metabolic_gene$geneSymbols

# Draw Venn plot
library(ggvenn)

a <- list(
  "Metastase_DEGS" = DEGS_NAME,
  "Metabolic_Genes" = meta_genes
)

ggvenn(
  a,
  c("Metastase_DEGS", "Metabolic_Genes"),
  show_elements = FALSE,
  show_percentage = TRUE,
  digits = 3,
  fill_color = c("#EF8A43", "#4865A9"),
  fill_alpha = 0.5,
  stroke_color = "black",
  stroke_alpha = 1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "#75A4C9",
  set_name_size = 5,
  text_color = "black",
  text_size = 5,
  label_sep = ","
)

# Save Venn plot as PDF
pdf("VENN PLOT of DEGS and metabolic genes.pdf")
# Example plot (replace with actual plot if needed)
plot(1:10, 1:10, main = "Example Plot")
dev.off()