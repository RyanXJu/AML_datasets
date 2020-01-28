# https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("TCGAbiolinks")
devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")

library(TCGAbiolinks)


# lanch error in Rstudio, but can run in terminal R
query <- GDCquery(project = "TARGET-AML",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts"
                  )

GDCdownload(query = query)