## GSE131418 dataset preprocessing
## Set groups
## PCA analysis

#########################################################################
## GEO data download and preprocessing

rm(list=ls())
library(BiocManager)
library(GEOquery)
options(stringsAsFactors = FALSE)
library(limma)
library(GEOquery) # Load GEOquery R Package

# Set working directory to ".../raw data/Functional DEGs/GSE131418"
# Upload gene expression matrix
# Upload the downloaded file "series_matrix.txt"
eSet = getGEO(filename = "GSE131418_series_matrix.txt.gz", getGPL = FALSE) # Open the locally downloaded GEO file

library(data.table)
# Upload from the file ".txt.gz"
data <- read.table("GSE131418_Consortium_prim_met_GE_probe_level.txt.gz", header = TRUE, sep = "\t")
class(data)
data[1:10,1:10]

data2 = read.table("GSE131418_MCC_prim_met_GE_probe_level.txt.gz", header = TRUE, sep = "\t")
# Combine two data.tables into one gene expression matrix
exp = cbind(data, data2)

# Normalize between arrays
boxplot(exp) # Check the intensity of expression between arrays
range(exp)
exp = normalizeBetweenArrays(exp) 
range(exp) 

# Assess the need to normalize the expression matrix
# Determine whether log transformation is required
qx <- as.numeric(quantile(exp, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE))
LogC <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  exp[which(exp <= 0)] <- NaN
  exp <- log2(exp + 1) 
}

range(exp)

exp = as.data.frame(exp) # Transform matrix to dataframe

# Acquire sample clinical information matrix 
ph = pData(eSet)
ph = as.data.frame(ph)

# -----------------------------------------------------------------------------

## Match Gene Symbol labels to gene IDs 
# Upload list of gene symbols from "GSE 131418 symbol.xlsx"
library(readxl)
ID = read_excel("GSE 131418 symbol.xlsx", sheet = 1, col_names = TRUE) 
ID = as.data.frame(ID)

# Remove rows containing NA in GeneSymbol
ID <- ID[!is.na(ID$GeneSymbol), ]

# Transfer IDs to Gene Symbol
exp$ID = rownames(exp)
exp = merge(exp, ID, by = "ID")
exp = as.data.frame(exp)
exp[1:3,1135:1140]
exp = exp[, -c(1137, 1138, 1140)]
exp[1:5,1:5]
exp = exp[, -1]
rownames(exp) = exp$GeneSymbol

# Remove duplicated gene symbols in the gene expression dataframe
which(duplicated(exp$GeneSymbol))
k = !duplicated(exp$GeneSymbol)
table(k)
exp = exp[k,]

# Copy gene symbol to row names of the dataframe
rownames(exp) = exp$GeneSymbol
colnames(exp)[1130:1136]
ll = head(exp)
exp = exp[,-1136]

##################
## Data filter
# Remove the samples with medical treatment
table(ph```r
      treatment status classification (see description):ch1`)
ph = ph[ph```r
        treatment status classification (see description):ch1` %in% "PRE", ]

# Match samples without medical treatment to gene expression dataframe
k = intersect(colnames(exp), ph``r
              treatment status classification (see description):ch1`)
ph = ph[ph```r
        treatment status classification (see description):ch1` %in% "PRE", ]

# Match samples without medical treatment to gene expression dataframe
k = intersect(colnames(exp), ph$title)
ph = ph[ph$title %in% k, ]
exp = exp[, k]
identical(ph$title, colnames(exp))

# Reassure sample order in gene expression dataframe
index = match(ph$title, colnames(exp))
exp = exp[, index]
identical(ph$title, colnames(exp))

## Save preprocessed data from GSE131418 as .Rdata
save(exp, ID, ph, file = "GSE131418_Raw_Data.Rdata")

###################################################################

## Set groups
# Labels : metastasis & Primary CRC
# Summary of labels 
table(ph$characteristics_ch1.6)
# Set group of samples to metastasis & primary
grouplist = ifelse(ph$characteristics_ch1.6 == "tumor type: METASTASIS", "metastasis", "primary")

# Set grouping labels as factors
grouplist = factor(grouplist, levels = c("metastasis", "primary"))
table(grouplist)

# Save sample groupings 
save(grouplist, file = "group_metastasis_primary.Rdata")

###################################################################
## PCA Analysis

rm(list=ls())

# Load data of samples and groupings
load(file = "group_Metastase.Rdata")
load(file = "GSE131418_Raw_Data.Rdata")

library(limma)
library(ggplot2)

# Data processing
dat = t(exp)
dat[is.na(dat)] = 0
df_pca = prcomp(dat)
df_pcs = data.frame(df_pca$x, Species = grouplist) 
head(df_pcs, 3)

# Visualization 
plot(df_pca$x[,1], df_pca$x[,2])
ggplot(df_pcs, aes(x = PC1, y = PC2, color = Species)) + geom_point()

# Modified PCA Plot
ggplot(df_pcs, aes(x = PC1, y = PC2, color = Species)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# Further PCA Plot customization
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste(colnames(df_pcs), "(", paste(as.character(percentage), "%", ")", sep = ""))
ggplot(df_pcs, aes(x = PC1, y = PC2, color = Species)) +
  geom_point() + 
  xlab(percentage[1]) +
  ylab(percentage[2])

# Final PCA Plot modifications
ggplot(df_pcs, aes(x = PC1, y = PC2, color = Species)) +
  geom_point() +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  theme(panel.background = element_blank(),  # Remove background color
        panel.grid = element_blank())        # Remove grid lines

## Save plot as PDF in the working directory
pdf("PCA Analysis.pdf")
# Generate illustration
plot(1:10, 1:10, main = "Example Plot")
# Close PDF device
dev.off()[ph```r
          treatment status classification (see description):ch1` %in% "PRE", ]

# Match samples without medical treatment to gene expression dataframe
k = intersect(colnames(exp), ph$title)
ph = ph[ph$title %in% k, ]
exp = exp[, k]
identical(ph$title, colnames(exp))

# Reassure sample order in gene expression dataframe
index = match(ph$title, colnames(exp))
exp = exp[, index]
identical(ph$title, colnames(exp))

## Save preprocessed data from GSE131418 as .Rdata
save(exp, ID, ph, file = "GSE131418_Raw_Data.Rdata")

###################################################################

## Set groups
# Labels : metastasis & Primary CRC
# Summary of labels 
table(ph$characteristics_ch1.6)
# Set group of samples to metastasis & primary
grouplist = ifelse(ph$characteristics_ch1.6 == "tumor type: METASTASIS", "metastasis", "primary")

# Set grouping labels as factors
grouplist = factor(grouplist, levels = c("metastasis", "primary"))
table(grouplist)

# Save sample groupings 
save(grouplist, file = "group_metastasis_primary.Rdata")

###################################################################
## PCA Analysis

rm(list=ls())

# Load data of samples and groupings
load(file = "group_Metastase.Rdata")
load(file = "GSE131418_Raw_Data.Rdata")

library(limma)
library(ggplot2)

# Data processing
dat = t(exp)
dat[is.na(dat)] = 0
df_pca = prcomp(dat)
df_pcs = data.frame(df_pca$x, Species = grouplist) 
head(df_pcs, 3)

# Visualization 
plot(df_pca$x[,1], df_pca$x[,2])
ggplot(df_pcs, aes(x = PC1, y = PC2, color = Species)) + geom_point()

# Modified PCA Plot
ggplot(df_pcs, aes(x = PC1, y = PC2, color = Species)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# Further PCA Plot customization
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste(colnames(df_pcs), "(", paste(as.character(percentage), "%", ")", sep = ""))
ggplot(df_pcs, aes(x = PC1, y = PC2, color = Species)) +
  geom_point() + 
  xlab(percentage[1]) +
  ylab(percentage[2])

# Final PCA Plot modifications
ggplot(df_pcs, aes(x = PC1, y = PC2, color = Species)) +
  geom_point() +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  theme(panel.background = element_blank(),  # Remove background color
        panel.grid = element_blank())        # Remove grid lines

## Save plot as PDF in the working directory
pdf("PCA Analysis.pdf")
# Generate illustration
plot(1:10, 1:10, main = "Example Plot")
# Close PDF device
dev.off()