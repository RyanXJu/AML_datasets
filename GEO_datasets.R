####  Load datasets from GEO ####

######## install package ##########################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GEOquery")
# BiocManager::install("Biobase")

library(GEOquery)
library(Biobase)

folder = setwd("/u/juxiao/GEO_data")
getwd()

# GSE49642 43 AML rna-seq (Lavallée VP et al. 2013)
# GSE6891 461+76 AML microarray (Verhaak RG et al. 2008)
# GSE30652  239 stem cell microarray 	(Nazor KL et al. 2012)
#    labels: unique(pData(test[[1]])[, 12])
# GSE42414 41 human cord blood cell surface expression ()
#     prob to genes: https://www.biostars.org/p/109248/
# GSE24759 38 hematopoietic cell surface marker expression 
# GSe12417 156+70 CN-AML cohorts
# 

test <- getGEO("GSE76009",GSEMatrix=TRUE) 
data<-exprs(test[[1]])  # get GE data
pheno <- pData(test[[1]]) # get clinical data

library("biomaRt")
# ensembl = useMart(biomart= "ensembl",dataset="hsapiens_gene_ensembl")
# affy_ensembl= c("affy_hg_u133_plus_2", "ensembl_gene_id")
# getBM(attributes= affy_ensembl, mart= ensembl, values = "*", uniqueRows=T)

## transfer prob id to gene_id
## inspired by https://www.biostars.org/p/76097/
genes22hugo <- function( v, srcType = "affy_hg_u133_plus_2" )
{
  ## Retrieve the EMSEMBL -> HUGO mapping
  ensembl = useMart(biomart= "ensembl",dataset="hsapiens_gene_ensembl")
  ID <- biomaRt::getBM( attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl )
  
  ## Make sure there was at least one mapping
  if( nrow(ID) < 1 ) top( "No IDs mapped successfully" )
  
  ## Drop empty duds
  j <- which( ID[,2] == "" )
  if( length(j) > 0 ) ID <- ID[-j,]
  stopifnot( all( ID[,1] %in% v ) )
  
  ID
}


genes <- genes22hugo( rownames(data)[1:5] )
