### Function transform probe information to genes ids ###


# @scrType: source type of the probes
# listAttributes(ensembl) to see all possible types
probe2gene <- function( v, srcType = NULL )
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