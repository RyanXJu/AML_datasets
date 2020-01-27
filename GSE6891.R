# GSE6891
# 
# data: 461+76 AML microarray
# platform: affy_hg_u133_plus_2
# values: Log2-transformed MAS5 processed intensity after applying a minimum threshold of 30
# Origine: Verhaak, Roel GW, et al.2009
# Used in: Wang, Mei, et al. 2018

library(GEOquery)
library(Biobase)
library(biomaRt)

folder = setwd("/u/juxiao/AML-RNAseq_Datasets/GSE6891")
getwd()
fnGE = "GSE6891.tsv"
fnY = "GSE6891_labels.tsv"


gse <- getGEO("GSE6891",GSEMatrix=TRUE) 
data<-exprs(gse[[1]])  # get GE data
pheno <- pData(gse[[1]]) # get clinical data

data[1:3,1:3]
#            GSM158711 GSM158712 GSM158713
# 1007_s_at  4.906891  7.008989  4.906891
# 1053_at    7.071462  5.514122  6.757557
# 117_at     5.181898  5.459432  6.006747

platform <- "affy_hg_u133_plus_2"
source("probe2gene.R")

probes <- rownames(data)
# probes[1:3]
# [1] "1007_s_at" "1053_at"   "117_at"   


# Map Ensembl IDs to HUGO
gene_id <- probe2gene( rownames(data), srcType = platform  )
head(gene_id)
# affy_hg_u133_plus_2 hgnc_symbol
# 1        1553551_s_at      MT-ND1
# 2        1553551_s_at       MT-TI
# 3        1553551_s_at       MT-TM
# 4        1553551_s_at      MT-ND2
# 5        1553538_s_at      MT-CO1
# 6          1553569_at      MT-CO1
# problem: one probe can map to several genes
# TODO!!!!!!

X <- data[gene_id[,1],]
rownames(X) <- gene_id[,2]
dim(X)

Y <- pheno

write.table(X, file = fnGE, sep = "\t", quote = FALSE, col.names = TRUE)
write.table(Y, file = fnY, sep = "\t", quote = FALSE, col.names = TRUE)
