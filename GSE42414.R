# GSE42414
# 
# data: 41 human cord blood cell surface expression
# rna_seq: microarray
# platform: illumina_humanht_12_v4
# Origine: Laurenti, Elisa, et al. 2013
# Used in: Ng, Stanley WK, et al.2016


library(GEOquery)
library(Biobase)
library(biomaRt)

folder = setwd("/u/juxiao/AML_datasets")
getwd()
fnGE = "/u/juxiao/geo_datasets/GSE42414/GSE42414.tsv"
fnY = "/u/juxiao/geo_datasets/GSE42414/GSE42414_labels.tsv"


gse <- getGEO("GSE42414",GSEMatrix=TRUE) 
data<-exprs(gse[[1]])  # get GE data
pheno <- pData(gse[[1]]) # get clinical data

data[1:3,1:3]
#              GSM1039578 GSM1039579 GSM1039580
# ILMN_1343291  14.715726  14.603213  14.585952
# ILMN_1651209   9.146223   9.256855   8.520572
# ILMN_1651228  13.713320  13.901335  13.753404

platform <- "illumina_humanht_12_v4"
source("probe2gene.R")

probes <- rownames(data)
probes[1:3]
# [1] "ILMN_1343291" "ILMN_1651209" "ILMN_1651228"


# Map Ensembl IDs to HUGO
gene_id <- probe2gene( rownames(data), srcType = platform )
head(gene_id)
# illumina_humanht_12_v4 hgnc_symbol
# 1           ILMN_1653529       TEX10
# 2           ILMN_1653463     SULT1C3
# 3           ILMN_1651850       RPS16
# 4           ILMN_1658081      FRMPD2
# 5           ILMN_1651236      GGTLC1
# 6           ILMN_1658123        DAZ1


X <- data[gene_id[,1],]
rownames(X) <- gene_id[,2]
dim(X)

#Y <- pheno[,c("cell type:ch1")]
Y <- pheno

write.table(X, file = fnGE, sep = "\t", quote = FALSE, col.names = TRUE)
write.table(Y, file = fnY, sep = "\t", quote = FALSE, col.names = TRUE)
