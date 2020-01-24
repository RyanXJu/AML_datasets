# GSE76009
# data: 138lsc+ & 89lsc- cell fractions from 78 AML
# rna_seq: microarray
# platform: illumina_humanht_12_v4
# labels: cd34/cd38, lsc+/lsc-
# Origine: Ng, Stanley WK, et al.2016 

library(GEOquery)
library(Biobase)
library(biomaRt)

folder = setwd("/u/juxiao/AML-RNAseq_Datasets/GSE76009")
getwd()
fnGE = "GSE76009.tsv"
fnY = "GSE76009_labels.tsv"


gse <- getGEO("GSE76009",GSEMatrix=TRUE) 
data<-exprs(gse[[1]])  # get GE data
pheno <- pData(gse[[1]]) # get clinical data

data[1:3,1:3]
# GSM1972167 GSM1972168 GSM1972169
# ILMN_1343291   15.53806  15.451766  15.519726
# ILMN_1343295   12.69164  12.296734  13.378752
# ILMN_1651199    6.29835   6.369729   6.425103

platform <- "illumina_humanht_12_v4"
source("probe2gene.R")

probes <- rownames(data)
# probes[1:3]
# [1] "ILMN_1343291" "ILMN_1343295" "ILMN_1651199"


# Map Ensembl IDs to HUGO
gene_id <- probe2gene( rownames(data), srcType = platform )
head(gene_id)
# illumina_humanht_12_v4 hgnc_symbol
# 1           ILMN_1653529       TEX10
# 2           ILMN_1653463     SULT1C3
# 3           ILMN_1651850       RPS16
# 4           ILMN_1652784       CDHR5
# 5           ILMN_1651236      GGTLC1
# 6           ILMN_1652758        AGRN


X <- data[gene_id[,1],]
rownames(X) <- gene_id[,2]
dim(X)

Y <- pheno[,c("cd34/cd38 expression:ch1","lsc content:ch1")]

write.table(X, file = fnGE, sep = "\t", quote = FALSE, col.names = TRUE)
write.table(Y, file = fnY, sep = "\t", quote = FALSE, col.names = TRUE)

pheno
