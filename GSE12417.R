# GSE12417
# 
# data: 163+79 CN-AML cohorts
# rna_seq: microarray
# platform; Affymetrix Human Genome U133A & U133B (163 training), 
#            U133Plus2 (79 test)
# values: Normalization was performed using the variance stabilizing 
#           normalization algorithm,and probe set expression values 
#           were calculated by the median polish method.
# Origine: Metzeler, Klaus H., et al 2008
# Used in: Ng, Stanley WK, et al.2016


library(GEOquery)
library(Biobase)
library(biomaRt)

folder = setwd("/u/juxiao/AML_datasets")
getwd()
fnGE = "/u/juxiao/geo_datasets/GSE12417/GSE12417.tsv"
fnY = "/u/juxiao/geo_datasets/GSE12417/GSE12417_labels.tsv"


gse <- getGEO("GSE12417",GSEMatrix=TRUE) 

########## 79 test data ############
### platform : U133 Plus 2.0 (Affymetrix Human Genome U133 Plus 2.0 Array)
data_79<-exprs(gse[[1]])  # get GE data
pheno_79 <- pData(gse[[1]]) # get clinical data
data_79[1:3,1:3]
# GSM316652 GSM316653 GSM316654
# 1007_s_at  6.583821  6.437739  6.023578
# 1053_at    7.595103  8.919206  8.599511
# 117_at     7.976983  6.583684  6.779236

########## 163 training data part1############
### platform: U133A (Affymetrix Human Genome U133A Array)
data_163_1<-exprs(gse[[2]])  # get GE data
pheno_163_1 <- pData(gse[[2]]) # get clinical data
data_163_1[1:3,1:3]
# GSM311598 GSM311599 GSM311600
# 1007_s_at  9.855042  9.647614  9.308163
# 1053_at    7.851224  8.101409  8.326415
# 117_at     8.039294  8.069666  8.205875


########## 163 training data part2 ############
### platform: U133B (Affymetrix Human Genome U133B Array)
data_163_2<-exprs(gse[[3]])  # get GE data
pheno_163_2 <- pData(gse[[3]]) # get clinical data
data_163_2[1:3,1:3]
# GSM311761 GSM311762 GSM311763
# 222385_x_at  8.112904  8.530817  8.725246
# 222386_s_at  8.573234 10.257322  9.702850
# 222387_s_at  8.035575  8.612441  8.309333

# combine the information of the 163 training data
colnames(data_163_2)<- colnames(data_163_1)
new_163 <- rbind(data_163_1, data_163_2)

# combine testing and training sets
# remove platform specific probes from the two datasets
data_79 <- data_79[rownames(data_79) %in% rownames(new_163),]
new_163 <- new_163[rownames(new_163) %in% rownames(data_79),]
# in total 44754 probe kept

data <- cbind(new_163, data_79)
pheno <- unlist(list(pheno_163_1$characteristics_ch1, 
                     pheno_79$characteristics_ch1))
names(pheno) <- colnames(data)

platform <- "affy_hg_u133_plus_2"
source("probe2gene.R")

probes <- rownames(data)
probes[1:3]
# [1] "1007_s_at" "1053_at"   "117_at"


# Map Ensembl IDs to HUGO
gene_id <- probe2gene( rownames(data), srcType = platform )
head(gene_id)
# affy_hg_u133_plus_2 hgnc_symbol
# 1           201023_at        TAF7
# 2           201212_at        LGMN
# 3         201001_s_at      UBE2V1
# 4         201002_s_at      UBE2V1
# 5         201065_s_at     GTF2IP1
# 6           201258_at       RPS16


X <- data[gene_id[,1],]
rownames(X) <- gene_id[,2]
dim(X)

Y <- pheno

write.table(X, file = fnGE, sep = "\t", quote = FALSE, col.names = TRUE)
write.table(Y, file = fnY, sep = "\t", quote = FALSE, col.names = TRUE)



X <- data[gene_id[,1],]
rownames(X) <- gene_id[,2]
dim(X)

#Y <- pheno[,c("cell type:ch1")]
Y<- pheno

write.table(X, file = fnGE, sep = "\t", quote = FALSE, col.names = TRUE)
write.table(Y, file = fnY, sep = "\t", quote = FALSE, col.names = TRUE)