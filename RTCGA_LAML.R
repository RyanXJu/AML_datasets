# install.packages('survminer')
# 
# ## source("https://bioconductor.org/biocLite.R")
# ## Error: With R version 3.5 or greater, install Bioconductor packages using BiocManager; 
# ## see https://bioconductor.org/install
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# #BiocManager::install(version = "3.10")
# BiocManager::install()
# 
# 
# ## biocLite("RTCGA.clinical") # data for examples
# ## biocLite is no longer needed for R 3.6
# 
# install.packages("ggthemes")  ## dependencies of "RTCGA.clinical"
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("RTCGA.clinical")

### ggsurvplot example
## example https://cran.r-project.org/web/packages/survminer/vignettes/Informative_Survival_Plots.html
library(survminer)
library(RTCGA.clinical)
library(survival)


# BiocManager::install("RTCGA.rnaseq")
library(RTCGA.rnaseq) # genes' expression

# BiocManager::install("RTCGA.FPKM")


######  LAML.clinical,  n = 200 #####
survivalTCGA(LAML.clinical,
             extract.cols = "admin.disease_code") -> LAML.survInfo

###### LAML.rnaseq , n= 173 ####
## other patients don't have rna-seq information
expressionsTCGA(LAML.rnaseq) -> LAML.rnaseq
