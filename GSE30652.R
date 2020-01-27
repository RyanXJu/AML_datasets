# GSE30652
# 
# data: 239 stem cell microarray
# rna_seq: microarray
# platform: illumina_humanht_12_v4
# values: The data were normalized by RSN normalization using LUMI package in R
# Origine: Nazor et al. 2012
# Used in: Malta et al. 2018


library(GEOquery)
library(Biobase)
library(biomaRt)

folder = setwd("/u/juxiao/AML_datasets")
getwd()
fnGE = "/u/juxiao/geo_datasets/GSE30652/GSE30652.tsv"
fnY = "/u/juxiao/geo_datasets/GSE30652/GSE30652_labels.tsv"


gse <- getGEO("GSE30652",GSEMatrix=TRUE) 
data<-exprs(gse[[1]])  # get GE data
pheno <- pData(gse[[1]]) # get clinical data

data[1:3,1:3]
#              GSM759965   GSM759966   GSM759967
# ILMN_1343291 27562.72400 34471.36194 33229.61742
# ILMN_1343295 14289.37516 12244.94815 11113.99598
# ILMN_1651209    53.36326    53.73044    70.48524

platform <- "illumina_humanht_12_v4"
source("probe2gene.R")

probes <- rownames(data)
probes[1:3]
# [1] "ILMN_1343291" "ILMN_1343295" "ILMN_1651209"


# Map Ensembl IDs to HUGO
gene_id <- probe2gene( rownames(data), srcType = platform )
head(gene_id)
# illumina_humanht_12_v4 hgnc_symbol
# 1           ILMN_1653529       TEX10
# 2           ILMN_1651850       RPS16
# 3           ILMN_1651236      GGTLC1
# 5           ILMN_1656743        MIA2
# 6           ILMN_1657317      POLR2J
# 7           ILMN_1651886     CWF19L1


X <- data[gene_id[,1],]
rownames(X) <- gene_id[,2]
dim(X)

#Y <- pheno[,c("cell type:ch1")]
Y<- pheno

write.table(X, file = fnGE, sep = "\t", quote = FALSE, col.names = TRUE)
write.table(Y, file = fnY, sep = "\t", quote = FALSE, col.names = TRUE)