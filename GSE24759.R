# GSE24759
# 
# data: 38 hematopoietic cell surface marker expression
# rna_seq: microarray
# platform: illumina_humanht_12_v4
# values: log2-transformed RMA data and normalized such that the mean of each gene is zero
# Origine: Novershtern, Noa, et al. 2011
# Used in: Ng, Stanley WK, et al.2016


library(GEOquery)
library(Biobase)
library(biomaRt)

folder = setwd("/u/juxiao/AML_datasets")
getwd()
fnGE = "/u/juxiao/geo_datasets/GSE24759/GSE24759.tsv"
fnY = "/u/juxiao/geo_datasets/GSE24759/GSE24759_labels.tsv"


gse <- getGEO("GSE24759",GSEMatrix=TRUE) 
data<-exprs(gse[[1]])  # get GE data
pheno <- pData(gse[[1]]) # get clinical data

data[1:3,1:3]
#             GSM609632 GSM609633 GSM609634
# 1007_s_at    6.3116    6.2723    5.8962
# 1053_at      6.5486    6.3471    6.4733
# 117_at       6.2804    6.6932    6.3526

platform <- "affy_hg_u133_plus_2"
source("probe2gene.R")

probes <- rownames(data)
probes[1:3]
# [1] "1007_s_at" "1053_at"   "117_at" 


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