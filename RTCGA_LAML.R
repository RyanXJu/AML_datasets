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
library(RTCGA)
installTCGA("RTCGA.clinical")

# BiocManager::install("RTCGA.rnaseq")
library(RTCGA.rnaseq) # genes' expression

# BiocManager::install("RTCGA.FPKM")


folder = setwd("/u/juxiao/AML_datasets")
getwd()
fnGE = "/u/juxiao/synapse_datasets/syn4557014_TCGA/LAML_TCGA.tsv"
fnY = "/u/juxiao/synapse_datasets/syn4557014_TCGA/LAML_TCGA_labels.tsv"
fd = "/u/juxiao/synapse_datasets/syn4557014_TCGA/"


###### LAML.rnaseq , n= 173 ####
## other patients don't have rna-seq information
LAML.rnaseq <- expressionsTCGA(LAML.rnaseq)


###### LAML.clinical, n=200 #######
keeps <- c("patient.bcr_patient_barcode",
           "patient.acute_myeloid_leukemia_calgb_cytogenetics_risk_category",
           "patient.age_at_initial_pathologic_diagnosis",
           "patient.leukemia_french_american_british_morphology_code",
           "patient.vital_status",
           "patient.days_to_last_followup",
           "patient.days_to_death"
           )

clinical <- LAML.clinical[, colnames(LAML.clinical) %in% keeps]

# ######  LAML.survInfo,  n = 200 #####
# ###### survival time = last follow up/ death
# LAML.survInfo <- survivalTCGA(LAML.clinical,
#                               extract.cols = "admin.disease_code")

X<- LAML.rnaseq
Y<- clinical

write.table(X, file = fnGE , sep = "\t", quote = FALSE, col.names = TRUE)
write.table(Y, file = fnY , sep = "\t", quote = FALSE, col.names = TRUE)
