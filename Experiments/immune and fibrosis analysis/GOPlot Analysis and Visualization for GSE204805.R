
############## GOPlot Analysis and Visualization for GSE204805 #################
################################################################################

# Section 1: Setup and Data Preparation

library(GOplot)
library(dplyr)

# Format GO enrichment results for GOplot
result_go_formatted <- result_go %>%
  mutate(
    Category = ONTOLOGY,
    ID = ID,
    Term = Description,
    Genes = gsub("/", ", ", geneID),
    adj_pval = p.adjust
  ) %>%
  select(Category, ID, Term, Genes, adj_pval)

# Prepare DEG data
nrDEG$ID <- rownames(nrDEG)
nrDEG <- nrDEG[, c(7, 1:6)]

# Create circle data for GOplot
circ <- circle_dat(result_go_formatted, nrDEG)
save(circ, file = "GOPlot circ of ACMSD high and low.rdata")

# Section 2: GOPlot Visualizations

# Bar plot for Cellular Component (CC) category
GOBar(subset(circ, category == 'CC'))

# Bubble plots
GOBubble(circ, labels = 20)
GOBubble(
  circ,
  title = 'Bubble plot with background colour',
  display = 'multiple',
  bg.col = TRUE,
  labels = 20
)

# Sort by adjusted p-value and scale count column
circ <- circ[order(circ$adj_pval), ]
circ$count <- circ$count / 5

# Bubble plot with additional options
GOBubble(
  circ,
  title = 'Bubble plot with background colour',
  display = 'multiple',
  bg.col = TRUE,
  ID = TRUE,
  table.col = TRUE,
  table.legend = FALSE,
  labels = 22
)

# Circle plots for top terms or specific GO IDs
GOCircle(circ, nsub = 10)
IDs <- c(
  'GO:0036037', 'GO:0050852', 'GO:0050853', 'GO:0042267', 'GO:0001909',
  'GO:0070942', 'GO:0070098', 'GO:0004896', 'GO:0032760', 'GO:0050729'
)
GOCircle(circ, nsub = IDs)

# Section 3: Chord and Heatmap Visualizations

# Prepare gene and process data for chord plot
genes <- circ[circ$term %in% c("extracellular matrix organization", "collagen-containing extracellular matrix"), ]
gene <- genes[, c("genes", "logFC")]
colnames(gene) <- c("ID", "logFC")
process <- c("extracellular matrix organization", "collagen-containing extracellular matrix")

# Create chord data
chord <- chord_dat(circ, gene, process)

# Filter for marker genes
markers <- c("CCR2", "CCR5", "CCR7", "CCL17", "CXCL9", "CD160", "CD28", "IDO1", "CTLA4", "CD274", "LAG3")
gene <- gene[gene$ID %in% markers, ]
gene <- gene[!duplicated(gene$ID), ]

# Chord plot
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)

# Heatmap plots
GOHeat(chord[, -8], nlfc = 0)
GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))

# Cluster plot
GOCluster(circ, process, clust.by = 'logFC', term.width = 2)
GOCluster(circ, process, clust.by = 'logFC', lfc.col = c("darkred", "white", "lightblue"), term.width = 2)