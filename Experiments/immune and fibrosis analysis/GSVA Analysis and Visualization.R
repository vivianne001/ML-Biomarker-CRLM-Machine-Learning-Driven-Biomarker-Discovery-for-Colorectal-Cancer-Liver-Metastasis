
################ GSVA Analysis and Visualization for GSE204805 ################
################################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(GSVA)
library(GSEABase)
library(pheatmap)
library(limma)
library(BiocParallel)

# Section 1: Data Preparation

# load("GSE204805_processed.Rdata")

ph <- ph[ph$`cell type:ch1` == "LM", ]
k <- intersect(rownames(ph), colnames(exp))
ph <- ph[k, ]
exp <- exp[, k]

exp <- log2(exp + 1)
exp <- as.data.frame(exp)

# Section 2: Download and Prepare Gene Sets

# KEGG gene sets
KEGG_df_all <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
KEGG_df <- dplyr::select(KEGG_df_all, gs_name, gs_exact_source, gene_symbol)
kegg_list <- split(KEGG_df$gene_symbol, KEGG_df$gs_name)

# GO gene sets
GO_df_all <- msigdbr(species = "Homo sapiens", category = "C5")
GO_df <- dplyr::select(GO_df_all, gs_name, gene_symbol, gs_exact_source, gene_symbol)
go_list <- split(GO_df$gene_symbol, GO_df$gs_name)

save(KEGG_df, KEGG_df_all, kegg_list, file = "msigdbr-KEGG-GO-pathways.Rdata")

# Section 3: GSVA Calculation

dat <- as.matrix(exp)

# GSVA for GO gene sets
geneset <- go_list
gsva_mat <- gsva(
  expr = dat,
  gset.idx.list = geneset,
  kcdf = "Poisson",
  verbose = TRUE,
  parallel.sz = parallel::detectCores()
)
write.csv(gsva_mat, "gsva_go_matrix.csv")

# GSVA for KEGG gene sets
geneset <- kegg_list
gsva_mat <- gsva(
  expr = dat,
  gset.idx.list = geneset,
  kcdf = "Poisson",
  verbose = TRUE,
  parallel.sz = parallel::detectCores()
)
write.csv(gsva_mat, "gsva_kegg_matrix.csv")

gsva_mat <- read.csv(file = "gsva_kegg_matrix.csv", header = TRUE)

# Section 4: Differential Analysis with LIMMA

rownames(gsva_mat) <- gsva_mat$X
gsva_mat <- gsva_mat[, -1]

identical(rownames(ph), colnames(gsva_mat))
identical(rownames(ph), colnames(exp))

dat <- as.data.frame(t(exp))
median_val <- median(dat$ACMSD)
dat$grouplist <- ifelse(dat$ACMSD > median_val, "high", "low")
table(dat$grouplist)
grouplist <- dat$grouplist

exp_group <- "high"
ctr_group <- "low"

design <- model.matrix(~0 + factor(grouplist))
design <- as.data.frame(design)
rownames(design) <- colnames(gsva_mat)
colnames(design) <- levels(factor(grouplist))

contrast.matrix <- makeContrasts(contrasts = paste0(exp_group, '-', ctr_group), levels = design)

fit1 <- lmFit(gsva_mat, design)
fit2 <- contrasts.fit(fit1, contrast.matrix)
efit <- eBayes(fit2)

summary(decideTests(efit, lfc = 1, p.value = 0.05))
tempOutput <- topTable(efit, coef = paste0(exp_group, '-', ctr_group), n = Inf)
degs <- na.omit(tempOutput)

save(degs, file = "DEG of gsva_go_matrix.Rdata")

# Section 5: Heatmap Visualization

degs$ID <- rownames(degs)
degs <- degs[grepl("KEGG_", degs$ID), ]

padj_cutoff <- 0.05
log2FC_cutoff <- 0.1

keep <- rownames(degs[degs$adj.P.Val < padj_cutoff & abs(degs$logFC) > log2FC_cutoff, ])
dat_heatmap <- degs[keep, ]
dat_heatmap <- dat_heatmap[order(abs(dat_heatmap$logFC), decreasing = TRUE), ]
datt <- gsva_mat[rownames(dat_heatmap), ]

design <- as.data.frame(sapply(design, as.numeric))

library(RColorBrewer)
ph$group <- grouplist
ann_col <- ph[, c("geo_accession", "group")]
ann_col <- as.data.frame(ann_col)
ann_col$group <- factor(ann_col$group, levels = c("high", "low"))
ann_col <- ann_col[order(ann_col$group), ]
cols <- rownames(ann_col)
ann_col <- ann_col[, -1]
rownames(ann_col) <- cols
colnames(ann_col) <- "Group"
ann_col$Group <- as.factor(ann_col$Group)
ann_color <- list(Group = c(high = "#0089CF", low = "#E889BD"))

k <- rownames(ann_col)
datt <- datt[, k]

rownames(datt) <- sub("KEGG_", "", rownames(datt))
write.csv(datt, file = "GSVA-GO-Metabolism.csv")

datt <- read.csv("GSVA-GO-Metabolism.csv")
rownames(datt) <- datt$X
datt <- datt[, -1]

new_order <- c("KEGG_TGF_BETA_SIGNALING_PATHWAY", "KEGG_ECM_RECEPTOR_INTERACTION")
datt <- datt[new_order, ]

ann_row <- data.frame(pathway = factor(c(rep("KEGG_TGF_BETA_SIGNALING_PATHWAY", 1), rep("KEGG_ECM_RECEPTOR_INTERACTION", 1))))
rownames(ann_row) <- rownames(datt)
unique_pathways <- unique(ann_row$pathway)
pathway_colors <- c(
  "KEGG_TGF_BETA_SIGNALING_PATHWAY" = "#AF8CBB",
  "KEGG_ECM_RECEPTOR_INTERACTION" = "#E08D8B"
)
ann_color <- list(pathway = pathway_colors)

pheatmap(
  datt, scale = "row",
  cluster_rows = FALSE, cluster_cols = FALSE,
  treeheight_col = 30, treeheight_row = 30,
  border_color = "grey60",
  fontsize = 10, fontsize_row = 6, fontsize_col = 10,
  show_colnames = FALSE, show_rownames = TRUE,
  main = "Heatmap of High & Low Expressional Levels of ACMSD in mCRC LM Samples",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  annotation_col = ann_col, annotation_row = ann_row,
  annotation_colors = ann_color,
  annotation_legend = TRUE,
  annotation_names_col = TRUE, annotation_names_row = TRUE
)

# Section 6: Bar Plot Visualization of T Values

load("DEG_of_gsva_go_matrix.Rdata")

dat_plot <- data.frame(id = rownames(degs), t = degs$t)
library(stringr)
dat_plot$id <- str_replace(dat_plot$id, "HP_", "")
dat_plot$id <- str_replace(dat_plot$id, "GOCC_", "")
dat_plot$id <- str_replace(dat_plot$id, "GOBP_", "")
dat_plot$id <- str_replace(dat_plot$id, "GOMF_", "")
dat_plot <- dat_plot[!duplicated(dat_plot$id), ]

dat_plot$threshold <- factor(
  ifelse(dat_plot$t > -5, ifelse(dat_plot$t >= 5, 'Up', 'NoSignifi'), 'Down'),
  levels = c('Up', 'Down', 'NoSignifi')
)
dat_plot <- dat_plot %>% arrange(t)
dat_plot$id <- factor(dat_plot$id, levels = dat_plot$id)

library(ggthemes)
library(ggprism)

p <- ggplot(data = dat_plot, aes(x = id, y = t, fill = threshold)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c('Up' = '#36638a', 'NoSignifi' = '#cccccc', 'Down' = '#7bcd7b')) +
  geom_hline(yintercept = c(-5, 5), color = 'white', size = 0.5, lty = 'dashed') +
  xlab('') +
  ylab('t value of GSVA score, High and Low Expressional Levels of ACMSD in mCRC LM') +
  guides(fill = FALSE) +
  theme_prism(border = TRUE) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p