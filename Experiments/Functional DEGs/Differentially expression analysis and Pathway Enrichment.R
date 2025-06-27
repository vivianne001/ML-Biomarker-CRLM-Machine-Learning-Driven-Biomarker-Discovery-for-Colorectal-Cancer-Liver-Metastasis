###################### Differential Expression and Pathway Enrichment Analysis of Functional DEGs #########################
# GSE131418
rm(list = ls())

## Differential Expression Analysis
#######################################################################

# Load raw data and grouping information
load(file = "GSE131418_Raw_Data.Rdata")
load(file = "group_metastasis_primary.Rdata")

# Load required libraries
options(stringsAsFactors = FALSE)
library(ggplot2)
library(ggsignif)
library(ggsci)
library(stringr)
library(limma)
library(data.table)
library(maftools)
library(tidyverse)

# Summary of grouping
table(grouplist)

# Set experimental and control groups
expt = "metastasis"
ctr = "primary"

# Design matrix for Limma DEG analysis
design = model.matrix(~0 + factor(grouplist))
design = as.data.frame(design)
colnames(design) = levels(factor(grouplist))
rownames(design) = ph$title
design = design[order(design$metastasis, decreasing = TRUE), ]

ll = match(rownames(design), colnames(exp))
exp = exp[, ll]
ph = ph[ll, ]

colnames(design) = levels(factor(grouplist))
rownames(design) = colnames(exp)

contrast.matrix = makeContrasts(contrasts = paste0(expt, '-', ctr), levels = design)

fit1 = lmFit(exp, design)
fit2 = contrasts.fit(fit1, contrast.matrix)
efit = eBayes(fit2)

# Summary of DEGs
summary(decideTests(efit, lfc = 1, p.value = 0.05))
tempOutput = topTable(efit, coef = paste0(expt, '-', ctr), n = Inf)
nrDEG = na.omit(tempOutput)

range(nrDEG$logFC)
which(is.na(nrDEG))

nrDEG$gene = rownames(nrDEG)
nrDEG$gene = sapply(strsplit(nrDEG$gene, " "), `[`, 1)

write.csv(nrDEG, "nrDEG_grouping_results.csv")

range(nrDEG$logFC)
range(nrDEG$P.Value)

################################################

## Volcano Plot
LogFC = 1
P.value = 0.05
Q.value = 0.05
k1 = (nrDEG$P.Value < 0.05) & (nrDEG$logFC > 1)
k2 = (nrDEG$P.Value < 0.05) & (nrDEG$logFC < -1)
change = ifelse(k1, "up", ifelse(k2, "down", "stable"))
table(change)
nrDEG = mutate(nrDEG, change)

plot(nrDEG$logFC, -log10(nrDEG$P.Value))
DEG = nrDEG

logFC_cutoff = 1
this_tile = paste0(
  'Cutoff for logFC is ', round(logFC_cutoff, 3),
  '\nThe number of up gene is ', nrow(DEG[DEG$change == 'up', ]),
  '\nThe number of down gene is ', nrow(DEG[DEG$change == 'down', ])
)
head(DEG)

g = ggplot(data = DEG, aes(x = logFC, y = -log10(P.Value), color = change)) +
  geom_point(alpha = 0.5, size = 1.5) +
  theme_set(theme_set(theme_bw(base_size = 20))) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle(this_tile) + theme(plot.title = element_text(size = 15, hjust = 0.5)) +
  scale_colour_manual(values = c('blue', 'black', 'red')) +
  theme_bw()

print(g)

g = ggplot(data = DEG, aes(x = logFC, y = -log10(P.Value), color = change)) +
  geom_point(alpha = 0.5, size = 1.5) +
  theme_set(theme_set(theme_bw(base_size = 20))) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle(this_tile) + theme(plot.title = element_text(size = 15, hjust = 0.5)) +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black") +
  geom_hline(yintercept = -log10(P.value), lty = 4, col = "black") +
  scale_colour_manual(values = c('blue', 'black', 'red')) +
  theme_bw()

print(g)

library(magrittr)
DEG$gene_name = rownames(DEG)
if (TRUE) {
  x1 = DEG %>% filter(DEG$change == "up") %>% head(20)
  x2 = DEG %>% filter(DEG$change == "down") %>% head(20)
  for_label = rbind(x1, x2)
  for_label$gene_name = rownames(for_label)
}
for_label = as.data.frame(for_label)

volcano_plot = g +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = gene_name),
    data = for_label,
    color = "red"
  ) +
  theme_bw()

volcano_plot

save(for_label, nrDEG, DEG, file = "grouplist-volcano-metastasis-DATA.P.Value 0.05.Rdata")
ggsave(plot = g, file = "volcano_plot-INTERSECTGENES.pdf", height = 10, width = 12)

#################################################################################
## Pathway Enrichment
library(clusterProfiler)
library(hugene10sttranscriptcluster.db)
ids = toTable(hugene10sttranscriptclusterSYMBOL)

library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
library(tidyverse)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

DEG = read.csv(file = "GSE131418-DEGS.csv")
deg = DEG

load("functional dataset.Rdata")

s2e = bitr(DEG$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
deg = inner_join(deg, s2e, by = c("gene" = "SYMBOL"))
DEG = deg

gene_up = DEG[DEG$change == "up", "ENTREZID"]
gene_down = DEG[DEG$change == "down", "ENTREZID"]
gene_diff = c(gene_up, gene_down)

################################################################################ 
## GO Enrichment Plotting
ego = enrichGO(gene = gene_diff, OrgDb = org.Hs.eg.db, ont = "ALL", readable = TRUE)

dotplot(ego, split = "ONTOLOGY", font.size = 12, showCategory = 10) +
  facet_grid(ONTOLOGY ~ ., scale = "free") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 55)) +
  scale_fill_gradient(low = "red", high = "blue") +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 8))

ggsave(filename = "GO Enrichment of functional DEGs.pdf", width = 8, height = 10)