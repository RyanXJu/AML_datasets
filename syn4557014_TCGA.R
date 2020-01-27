



library(gelnet)
library(dplyr)
library(biomaRt)
library(synapser)

folder = setwd("/u/juxiao/AML_datasets")
getwd()
fnGE = "/u/juxiao/synapse_datasets/syn4557014_TCGA/syn4557014_TCGA.tsv"
fnY = "/u/juxiao/synapse_datasets/syn4557014_TCGA/syn4557014_TCGA_labels.tsv"
fd = "/u/juxiao/synapse_datasets/syn4557014_TCGA/"

# login to https://www.synapse.org
synLogin(email="xiaoju", password="02160505")


## load in pancan data as test
s <- synGet( "syn4976369", downloadLocation = getwd())

## Reduces HUGO|POSITION gene IDs to just HUGO
f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )

## preprocessing data, ** this may take 10-20minutes **
pancan <- read.delim( "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv", 
                      as.is=TRUE, check.names=FALSE ) %>%  ## Read the raw values
  filter( !grepl( "\\?", gene_id ) ) %>%      ## Drop genes with no mapping to HUGO
  mutate( gene_id = f( gene_id ) )      ## Clip gene ids to HUGO

write.table(pancan, file = fnGE , sep = "\t", quote = FALSE, col.names = TRUE)
