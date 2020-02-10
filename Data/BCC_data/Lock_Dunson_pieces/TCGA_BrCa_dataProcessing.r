###Adapted from the data processing file accompanying:
###Lock, E. F., & Dunson, D. B. (2013). Bayesian consensus clustering. Bioinformatics (Oxford, England), 29(20), 2610–2616. http://doi.org/10.1093/bioinformatics/btt425


# setwd("~/Dropbox/GwenDataIntegration/BCC")
# source("https://bioconductor.org/biocLite.R")

# if (!requireNamespace("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")

# BiocManager::install("impute")

library(impute)
library("MCMCpack")
library(magrittr)
library(stringr)

###Data can be downloaded from the TCGA Data Portal at: https://tcga-data.nci.nih.gov/docs/publications/brca_2012/

setwd("~/Documents/PhD/Year_1/Consensus_clustering/Data/BCC_data/")
GE      <- read.csv("BRCA.exp.348.med.csv", header = TRUE)
miRNA   <- read.csv("BRCA.348.precursor.txt", header = TRUE)
Protein <- read.table("rppaData-403Samp-171Ab-Trimmed.txt", header = TRUE)
Meth    <- read.table("BRCA.Methylation.574probes.802.txt", header = TRUE)

################Match columns (samples) between sources##############################
namesExp     <- names(GE)[2:349]
namesmiRNA   <- names(miRNA)[2:349]
namesProtein <- names(Protein)[2:404]
namesMeth    <- names(Meth)

# Matching samples present
namesExp     <- substr(namesExp,1,16)
namesmiRNA   <- substr(namesmiRNA,1,16)
namesProtein <- substr(namesProtein,1,16)

MatchProt    <- match(namesExp,namesProtein,nomatch=0)
MatchMeth    <- match(namesExp,namesMeth   ,nomatch=0)

miRNA_names  <- miRNA[,1]
miRNA        <- miRNA[,2:349]
miRNA.mat    <- as.matrix(miRNA[,order(namesmiRNA)]) %>% 
  set_rownames(miRNA_names)

Protein.mat  <- Protein[,2:404]
Protein.mat  <- as.matrix(Protein.mat[,MatchProt]) %>% 
  set_rownames(Protein[,1])

Meth.mat     <- as.matrix(Meth[,MatchMeth]) %>% 
  set_rownames(row.names(Meth))

###################Data processing#############################
set.seed(1)
###Impute missing expression values
#load package 'impute' from CRAN
Exp.mat      <- as.matrix(GE[,2:349])
Exp.mat      <- impute.knn(Exp.mat) ##Impute missing values via KNN (K=10)
Exp.mat      <- Exp.mat$data %>% 
  set_rownames(GE[,1])

rowSums(miRNA.mat==0)
##Remove miRNAs with > 50% missing values
miRNA.mat    <- miRNA.mat[rowSums(miRNA.mat==0) < 348*0.5,]

processedExpression  <- Exp.mat[apply(Exp.mat,1,sd)>1.5,] ###Filter to select only most variable genes
processedMethylation <- sqrt(Meth.mat)    ##take square root of methylation data
processedmiRNA       <- log(1+miRNA.mat) ##take log of miRNA data
processedProtein     <- scale(Protein.mat,center=TRUE,scale=TRUE) #Column center/scale protein

# Ask Paul re: this step
MDI_processed_Expression   <- processedExpression %>% 
  scale() %>% 
  t()

MDI_processed_Methylation  <- processedMethylation %>% 
  scale() %>% 
  t()

MDI_processed_miRNA        <- processedmiRNA %>% 
  scale() %>% 
  t()

MDI_processed_Protein      <- processedProtein %>% 
  scale() %>% 
  t()

# Participants should share row indices across datasets.
# Based on https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/ we want
# the participant component of the TGCA code. 
# Codes are of the format "TCGA.E2.A1B6.01A.21.A13E.20"; we want the "A1B6" from 
# this example.
# Use a REGEX pattern to extract this.
pattern <- "(?<=.)\\b(\\d|\\w){4}\\b"

# Check the participant order in each dataset.
participant_list <- vector("list", 4) %>% 
  set_names(c("Expression", "Protein", "Methylation", "miRNA"))

participant_list$Expression <- stringr::str_extract(row.names(MDI_processed_Expression), pattern)
participant_list$Protein <- stringr::str_extract(row.names(MDI_processed_Protein), pattern)
participant_list$miRNA <- stringr::str_extract(row.names(MDI_processed_miRNA), pattern)
participant_list$Methylation<- stringr::str_extract(row.names(MDI_processed_Methylation), pattern)

n <- nrow(MDI_processed_Expression)

# Do all pairwise combinations (unnecessary, if dataset 1 matches all others 
# they all match)
for(i in 1:3){
  for(j in (i+1):4){
    comp <- match(participant_list[[i]], participant_list[[j]])
    if(any(comp != 1:n)){
      cat("\nMismatch in order for datasets:\n")
      cat(names(participant_list)[i])
      cat(names(participant_list)[j])
      stop("Error.")
    }
  }
}

# Use the shortened names
# any(duplicated(participant_list$Expression))
# any(duplicated(participant_list$Protein))
# any(duplicated(participant_list$Methylation))
# any(duplicated(participant_list$miRNA))

row.names(MDI_processed_Expression) <- participant_list$Expression
row.names(MDI_processed_Methylation) <- participant_list$Methylation
row.names(MDI_processed_miRNA) <- participant_list$miRNA
row.names(MDI_processed_Protein) <- participant_list$Protein

# Save the MDI ready to be used by MDI
write.csv(MDI_processed_Expression, "./MDI_data/MDI_Expression_data.csv")
write.csv(MDI_processed_Methylation, "./MDI_data/MDI_Methylation_data.csv")
write.csv(MDI_processed_miRNA, "./MDI_data/MDI_miRNA_data.csv")
write.csv(MDI_processed_Protein, "./MDI_data/MDI_Protein_data.csv")

save(processedExpression , file = 'processedExpression.RDa')
save(processedMethylation, file = 'processedMethylation.RDa')
save(processedmiRNA      , file = 'processedmiRNA.RDa')
save(processedProtein    , file = 'processedProtein.RDa')
