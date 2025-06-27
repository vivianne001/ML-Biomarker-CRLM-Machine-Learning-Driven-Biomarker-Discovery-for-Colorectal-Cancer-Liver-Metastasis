
######################## Bulk GSEA Analysis Using KEGG Database ########################
# GSE204805
# Clear workspace
rm(list = ls())

# Load DEG results
load("ACMSD group_high and low.csv")
DEG = nrDEG

# Prepare gene list for GSEA
geneList = DEG$logFC
geneList = as.numeric(geneList)
names(geneList) = DEG$ENTREZID
geneList = sort(geneList, decreasing = TRUE)

# Remove duplicated SYMBOL values
DEG$symbol = rownames(DEG)
DEG$gene = DEG$symbol
DEG = DEG[!duplicated(DEG$gene), ]
rownames(DEG) = DEG$gene

# Load required libraries
library(clusterProfiler)
library(hugene10sttranscriptcluster.db)
ids = toTable(hugene10sttranscriptclusterSYMBOL)

library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)

# Check package version
packageVersion('pheatmap')

# Remove decimal points from ENSG IDs if present
library(stringr)
rownames(DEG) <- str_sub(rownames(DEG), start = 1, end = 15)
DEG$symbol <- DEG$gene

# Map SYMBOL to ENSEMBL and ENTREZID
df <- bitr(DEG$symbol,
           fromType = "SYMBOL",
           toType = c("ENSEMBL", "ENTREZID"),
           OrgDb = org.Hs.eg.db)

# Merge ENTREZID into DEG dataframe
colnames(DEG)[7] = "SYMBOL"
DEG <- inner_join(DEG, df, by = "SYMBOL")

# Prepare gene list for GSEA
gene_all <- as.character(DEG[, "ENTREZID"])
geneList <- DEG$logFC
names(geneList) <- DEG$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[geneList != 0]

# Run GSEA using KEGG
kk_gse <- gseKEGG(
  geneList     = geneList,
  organism     = 'hsa',
  minGSSize    = 5,
  pvalueCutoff = 0.95,
  verbose      = FALSE
)
res = kk_gse@result

save(res, kk_gse, file = "KEGG_Bulk_GSEA.Rdata")

############################## Visualization ##############################

# 1. GSEA Enrichment Plot for Selected Pathways
library(enrichplot)
library(ggplot2)

# Select target pathways
index = c("hsa04350", "hsa04064")
selected_gene_sets <- kk_gse@result$ID %in% index

# Plot GSEA enrichment curves
gseaplot2(
  kk_gse,
  geneSetID = kk_gse@result$ID[selected_gene_sets],
  ES_geom = 'line',
  pvalue_table = TRUE
)

# GSEA plot for a single pathway
g = gseaplot2(
  kk_gse,
  geneSetID = "hsa04512",
  color = "red",
  pvalue_table = TRUE,
  title = "ECM-receptor interaction"
)
g

############################## GO GSEA Analysis ##############################

# Run GSEA for GO Biological Process
egmt = gseGO(
  geneList = geneList,
  ont = "BP",
  OrgDb = "org.Hs.eg.db"
)
egmtd = egmt@result
egmt = na.omit(egmt)

save(egmtd, egmt, file = "GO_Bulk_GSEA.Rdata")

# Filter significant results
egmtd = egmtd[egmtd$NES > abs(1), ]
egmtd = egmtd[egmtd$p.adjust < 0.05, ]

# GSEA plot for a single GO term
g = gseaplot2(
  egmt,
  geneSetID = "GO:0032640",
  color = "red",
  pvalue_table = TRUE,
  title = "tumor necrosis factor production"
)
g

# GSEA plots for selected GO terms
index = c("GO:1990266", "GO:0002446")
selected_gene_sets <- egmt@result$ID %in% index
gseaplot2(
  egmt,
  geneSetID = egmt@result$ID[selected_gene_sets],
  ES_geom = 'line',
  pvalue_table = TRUE
)

# Visualize the 5th pathway
gseaplot2(
  egmt,
  5,
  color = "red",
  pvalue_table = TRUE,
  title = "DNA replication",
  base_size = 10,
  ES_geom = "line"
)

# Visualize multiple pathways in one plot
gseaplot2(
  egmt,
  1:3,
  color = "red",
  pvalue_table = TRUE,
  title = egmt@result$Description[1],
  base_size = 10,
  ES_geom = "line"
)

############################## Ridge Plot ##############################

# Ridge plot for core enrichment genes
selected_gene_sets <- kk_gse@result$ID %in% index
selected_gsea_result <- kk_gse@result[selected_gene_sets, ]
selected_gsea_result_obj <- new(
  "gseaResult",
  result = selected_gsea_result,
  geneSets = kk_gse@geneSets,
  geneList = kk_gse@geneList
)

ridgeplot(
  selected_gsea_result_obj,
  showCategory = 15,
  fill = "p.adjust",
  core_enrichment = TRUE,
  label_format = 30,
  orderBy = "NES",
  decreasing = FALSE
) +
  scale_fill_viridis_c()

############################## Pathway Network Plot ##############################

selected_gsea_result_obj = pairwise_termsim(selected_gsea_result_obj)
p1 = emapplot(selected_gsea_result_obj)
p2 = emapplot(selected_gsea_result_obj, cex_category = 1.5)
p3 = emapplot(selected_gsea_result_obj, layout = "kk")
p4 = emapplot(selected_gsea_result_obj, cex_category = 1.5, layout = "kk")

cowplot::plot_grid(p4, ncol = 1, labels = LETTERS[1:4])