# Load required packages
library(BiocManager)
library(GEOquery)
library(limma)
library(readxl)
library(FactoMineR)
library(factoextra)
library(data.table)
library(readr)

options(stringsAsFactors = FALSE)

############################## 1. GSE103479 ##############################
gse1 <- "GSE103479"
eSet1 <- getGEO(gse1, destdir=".", getGPL=FALSE)
exp1 <- exprs(eSet1[[1]])
exp1 <- as.matrix(exp1)
boxplot(exp1)
exp1 <- normalizeBetweenArrays(exp1)
# Log-transformation if necessary
qx <- as.numeric(quantile(exp1, c(0, .25, .5, .75, .99, 1.0), na.rm=TRUE))
logC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (logC) {
  exp1[exp1 <= 0] <- NA
  exp1 <- log2(exp1 + 1)
}
exp1 <- as.data.frame(exp1)
ph1 <- pData(eSet1[[1]])
ph1 <- as.data.frame(ph1)
ID1 <- read_excel("GPL23985-ID.xlsx", sheet = 1, col_names = TRUE)
ID1 <- as.data.frame(ID1)
ID1$`Gene Symbol` <- sub(" /// .*", "", ID1$`Gene Symbol`)
ID1$`Gene Symbol` <- make.unique(ID1$`Gene Symbol`)
ID1 <- ID1[!is.na(ID1$`Gene Symbol`), ]
exp1$ID <- rownames(exp1)
exp1 <- merge(exp1, ID1, by="ID")
rownames(exp1) <- exp1$`Gene Symbol`
exp1 <- exp1[,-c(ncol(exp1)-1, ncol(exp1))]
colnames(exp1)[110:112]
# Harmonize sample names
common1 <- intersect(colnames(exp1), rownames(ph1))
ph1 <- ph1[common1,]
exp1 <- exp1[, common1]
# Grouping
grouplist1 <- ifelse(ph1$`metastasis:ch1` == "none", "primary", "metastasis")
grouplist1 <- factor(grouplist1, levels = c("metastasis", "primary"))

save(exp1,ph1,file = "GSE103479_Raw_Data.Rdata")

############################## 2. GSE17536 ##############################
eSet2 <- getGEO(filename = "GSE17536_series_matrix.txt.gz", getGPL = FALSE)
exp2 <- exprs(eSet2)
exp2 <- as.matrix(exp2)
exp2 <- normalizeBetweenArrays(exp2)
qx <- as.numeric(quantile(exp2, c(0, .25, .5, .75, .99, 1.0), na.rm=TRUE))
logC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (logC) {
  exp2[exp2 <= 0] <- NA
  exp2 <- log2(exp2 + 1)
}
exp2 <- as.data.frame(exp2)
ph2 <- pData(eSet2)
ph2 <- as.data.frame(ph2)
ID2 <- read_excel("GSE39582 ID.xlsx", sheet = 1, col_names = TRUE)
ID2 <- as.data.frame(ID2)
ID2$`Gene Symbol` <- sub(" /// .*", "", ID2$`Gene Symbol`)
ID2$`Gene Symbol` <- make.unique(ID2$`Gene Symbol`)
ID2 <- ID2[!is.na(ID2$`Gene Symbol`), ]
exp2$ID <- rownames(exp2)
exp2 <- merge(exp2, ID2, by="ID")
rownames(exp2) <- exp2$`Gene Symbol`
exp2 <- exp2[,-c(ncol(exp2)-1, ncol(exp2))]
common2 <- intersect(colnames(exp2), rownames(ph2))
ph2 <- ph2[common2,]
exp2 <- exp2[, common2]
grouplist2 <- ifelse(ph2$`metastasis:ch1` == "none", "primary", "metastasis")
grouplist2 <- factor(grouplist2, levels = c("metastasis", "primary"))
save(exp2, ph2, file = "GSE17536_processed.Rdata")

############################## 3. GSE28702 ##############################
eSet3 <- getGEO(filename = "GSE28702_series_matrix.txt.gz", getGPL = FALSE)
exp3 <- exprs(eSet3)
exp3 <- as.matrix(exp3)
exp3 <- normalizeBetweenArrays(exp3)
qx <- as.numeric(quantile(exp3, c(0, .25, .5, .75, .99, 1.0), na.rm=TRUE))
logC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (logC) {
  exp3[exp3 <= 0] <- NA
  exp3 <- log2(exp3 + 1)
}
exp3 <- as.data.frame(exp3)
ph3 <- pData(eSet3)
ph3 <- as.data.frame(ph3)
ID3 <- read_excel("GSE39582 ID.xlsx", sheet = 1, col_names = TRUE)
ID3 <- as.data.frame(ID3)
ID3$`Gene Symbol` <- sub(" /// .*", "", ID3$`Gene Symbol`)
ID3$`Gene Symbol` <- make.unique(ID3$`Gene Symbol`)
ID3 <- ID3[!is.na(ID3$`Gene Symbol`), ]
exp3$ID <- rownames(exp3)
exp3 <- merge(exp3, ID3, by="ID")
rownames(exp3) <- exp3$`Gene Symbol`
exp3 <- exp3[,-c(ncol(exp3)-1, ncol(exp3))]
common3 <- intersect(colnames(exp3), rownames(ph3))
ph3 <- ph3[common3,]
exp3 <- exp3[, common3]
grouplist3 <- ifelse(ph3$`metastasis:ch1` == "none", "primary", "metastasis")
grouplist3 <- factor(grouplist3, levels = c("metastasis", "primary"))
save(exp3, ph3, file = "GSE28702_processed.Rdata")

############################## 4. GSE72970 ##############################
eSet4 <- getGEO(filename = "GSE72970_series_matrix.txt.gz", getGPL = FALSE)
exp4 <- exprs(eSet4)
exp4 <- as.matrix(exp4)
exp4 <- normalizeBetweenArrays(exp4)
qx <- as.numeric(quantile(exp4, c(0, .25, .5, .75, .99, 1.0), na.rm=TRUE))
logC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (logC) {
  exp4[exp4 <= 0] <- NA
  exp4 <- log2(exp4 + 1)
}
exp4 <- as.data.frame(exp4)
ph4 <- pData(eSet4)
ph4 <- as.data.frame(ph4)
ID4 <- read_excel("GSE39582 ID.xlsx", sheet = 1, col_names = TRUE)
ID4 <- as.data.frame(ID4)
ID4$`Gene Symbol` <- sub(" /// .*", "", ID4$`Gene Symbol`)
ID4$`Gene Symbol` <- make.unique(ID4$`Gene Symbol`)
ID4 <- ID4[!is.na(ID4$`Gene Symbol`), ]
exp4$ID <- rownames(exp4)
exp4 <- merge(exp4, ID4, by="ID")
rownames(exp4) <- exp4$`Gene Symbol`
exp4 <- exp4[,-c(ncol(exp4)-1, ncol(exp4))]
common4 <- intersect(colnames(exp4), rownames(ph4))
ph4 <- ph4[common4,]
exp4 <- exp4[, common4]
grouplist4 <- ifelse(ph4$`metastasis:ch1` == "none", "primary", "metastasis")
grouplist4 <- factor(grouplist4, levels = c("metastasis", "primary"))
save(exp4, ph4, file = "GSE72970_processed.Rdata")

############################## 5. GSE50760 ##############################
eSet5 <- getGEO(filename = "GSE50760_series_matrix.txt.gz", getGPL = FALSE)
exp5 <- read_excel("GSE50760-EXP.xlsx", col_names=T)
exp5 <- exp5[!duplicated(exp5$genes), ]
rownames(exp5) <- exp5$genes
exp5 <- exp5[,-1]
ph5 <- pData(eSet5)
ph5 <- as.data.frame(ph5)
ph5 <- ph5[ph5$source_name_ch1 %in% c("metastasized cancer", "primary colorectal cancer"), ]
common5 <- intersect(colnames(exp5), rownames(ph5))
exp5 <- exp5[, common5]
ph5 <- ph5[common5, ]
grouplist5 <- ifelse(ph5$source_name_ch1 == "metastasized cancer", "metastasis", "primary")
grouplist5 <- factor(grouplist5, levels = c("metastasis", "primary"))
save(exp5, ph5, file = "GSE50760_processed.Rdata")

############################## 6. GSE204805 ##############################
eSet6 <- getGEO(filename = "GSE204805_series_matrix.txt.gz", getGPL = FALSE)
data6 <- read_tsv("GSE204805_merged_hs_mm.tsv.gz")
data6 <- as.data.frame(data6)
data6$Geneid <- sub("H_", "", data6$Geneid)
data6$Geneid <- sapply(strsplit(data6$Geneid, "\\."), "[", 1)
data6 <- data6[!duplicated(data6$Geneid), ]
rownames(data6) <- data6$Geneid
exp6 <- data6[ , -1]
ph6 <- pData(eSet6)
ph6 <- as.data.frame(ph6)
ph6$title <- sapply(strsplit(ph6$title, ","), "[", 1)
rownames(ph6) <- ph6$title
exp6 <- exp6[, rownames(ph6)]
grouplist6 <- ifelse(ph6$`cell type:ch1` == "LM", "metastasis", "primary")
grouplist6 <- factor(grouplist6, levels = c("metastasis", "primary"))
save(exp6, ph6, grouplist6, file = "GSE204805_processed.Rdata")
save(grouplist6,file = "group_metastasis_primary.Rdata")

############################## 7. GSE108277 ##############################
eSet7 <- getGEO(filename = "GSE108277_series_matrix.txt.gz", getGPL = FALSE)
exp7 <- exprs(eSet7)
exp7 <- as.data.frame(exp7)
ph7 <- pData(eSet7)
ph7 <- as.data.frame(ph7)
ID7 <- read_excel("ID.xlsx", col_names = TRUE)
ID7 <- as.data.frame(ID7)
ID7 <- ID7[-c(1:13), ]
colnames(ID7) <- ID7[1, ]
ID7 <- ID7[-1, ]
ID7 <- ID7[complete.cases(ID7$ILMN_Gene), ]
ID7 <- ID7[!duplicated(ID7$ILMN_Gene), ]
exp7$ID <- rownames(exp7)
exp7 <- merge(exp7, ID7, by="ID")
index7 <- exp7$ILMN_Gene
exp7 <- exp7[,c(2:121)]
rownames(exp7) <- index7
# Normalization using the lumi package
if (!require("lumi")) BiocManager::install("lumi")
library(lumi)
exp7_matrix <- as.matrix(exp7)
eset7 <- new("ExpressionSet", exprs = exp7_matrix)
eset7_vst <- lumiT(eset7, method="log2")
eset7_norm <- lumiN(eset7_vst, method="quantile")
exp7_norm <- exprs(eset7_norm)
exp7 <- as.data.frame(exp7_norm)
save(exp7, ph7, file="GSE108277_processed.Rdata")

############################## 8. GSE159216 (example for annotation/normalization) ##############################
load("GSE159216_Raw_Data.Rdata")
exp8 <- as.data.frame(exp)
ID8 <- read_excel("GSE159216 ID.xlsx", sheet = 2, col_names=TRUE)
ID8 <- ID8[-c(1:14), ]
colnames(ID8) <- ID8[1, ]
ID8 <- ID8[-1, ]
ID8 <- ID8[, c(1,8)]
split_strings <- strsplit(ID8$gene_assignment, " // ")
ID8$gene_assignment <- sapply(split_strings, function(x) x[2])
bad_rows <- grepl("^LOC", ID8$gene_assignment) | ID8$gene_assignment == "---"
ID8 <- ID8[!bad_rows, ]
ID8 <- ID8[!is.na(ID8$gene_assignment), ]
exp8$ID <- rownames(exp8)
exp8 <- merge(exp8, ID8, by="ID")
exp8$`Gene Symbol` <- sub(" /// .*", "", exp8$`Gene Symbol`)
exp8 <- exp8[!duplicated(exp8$gene_assignment), ]
rownames(exp8) <- exp8$gene_assignment
exp8 <- exp8[,-c(ncol(exp8), 1)]
ph8 <- ph[rownames(exp8), ]
grouplist8 <- ifelse(ph8$`cell type:ch1` == "LM", "metastasis", "primary")
grouplist8 <- factor(grouplist8, levels = c("metastasis", "primary"))
save(exp8, ph8, file="GSE159216_processed.Rdata")