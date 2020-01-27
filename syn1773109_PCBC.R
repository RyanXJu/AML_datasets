# 
# ### Progenitor Cell Biology Consortium (PCBC) stem cell collection
# https://www.synapse.org/#!Synapse:syn1773109/wiki/54962
#   
#   Gene expression data: Normalized mRNA Matrix
# https://www.synapse.org/#!Synapse:syn2701943
#   
#   Stemness label data: RNA-Seq Metadata (301 samples, only 229 have gene expression data in the data above)
# https://www.synapse.org/#!Synapse:syn3156503/tables/

# values: Log2 Transform data and standardize



library(gelnet)
library(dplyr)
library(biomaRt)
library(synapser)

folder = setwd("/u/juxiao/AML_datasets")
getwd()
fnGE = "/u/juxiao/synapse_datasets/syn1773109_PCBC/syn1773109_PCBC.tsv"
fnY = "/u/juxiao/synapse_datasets/syn1773109_PCBC/syn1773109_PCBC_labels.tsv"
fd = "/u/juxiao/synapse_datasets/syn1773109_PCBC/"

# transfer probe info to gene_id
genes2hugo <- function( v, srcType = "ensembl_gene_id" )
{
  ## Retrieve the EMSEMBL -> HUGO mapping
  ensembl <- biomaRt::useMart( "ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="hsapiens_gene_ensembl" )
  ID <- biomaRt::getBM( attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl )
  
  ## Make sure there was at least one mapping
  if( nrow(ID) < 1 ) top( "No IDs mapped successfully" )
  
  ## Drop empty duds
  j <- which( ID[,2] == "" )
  if( length(j) > 0 ) ID <- ID[-j,]
  stopifnot( all( ID[,1] %in% v ) )
  
  ID
}

# login to https://www.synapse.org
synLogin(email="xiaoju", password="02160505")


###### Gene Expression data ########
# download GE data
synRNA <- synGet( "syn2701943", downloadLocation = fd )

X <- read.delim( paste(fd, "rnaseq_norm.tsv", sep = "") ) %>% 
  tibble::column_to_rownames( "tracking_id" ) %>% 
  as.matrix

X[1:3,1:3]
#                     H9.102.2.5 H9.102.2.6 H9.119.3.7
# ENSG00000000003.10  0.6306521  0.6071539  0.7197784
# ENSG00000000005.5  -1.2838794 -1.0152271 -0.8893850
# ENSG00000000419.8   0.6509436  0.5833117  0.7618753

## Drop the splice form ID from the gene names
v <- strsplit( rownames(X), "\\." ) %>% lapply( "[[", 1 ) %>% unlist()
rownames(X) <- v

# Map Ensembl IDs to HUGO
V <- genes2hugo( rownames(X) )
head(V)

X <- X[V[,1],]
rownames(X) <- V[,2]

#Reduce gene set to the provides list.
## fnGenes - [optional] filename of the list of entrez ID to consider
fnGenes = NULL
if(!is.null(fnGenes)){
  vGenes <- read.delim( fnGenes, header=FALSE ) %>% as.matrix() %>% drop()
  VE <- genes2hugo( vGenes, "entrezgene" )
  X <- X[intersect( rownames(X), VE[,2] ),]
}

cat(paste("Total number of genes in signuature: ", nrow(X), "\n"))

# mean centerized
# m <- apply( X, 1, mean )
# X <- X - m

########### label data #################
# download label data
synMeta <- synTableQuery( "SELECT UID, Diffname_short FROM syn3156503")
synMeta <- as.data.frame(synMeta)
Y<- synMeta

write.table(X, file = fnGE, sep = "\t", quote = FALSE, col.names = TRUE)
write.table(synMeta, file = fnY , sep = "\t", 
            quote = FALSE, col.names = TRUE)