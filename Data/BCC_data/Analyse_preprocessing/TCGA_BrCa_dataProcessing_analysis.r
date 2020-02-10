###Adapted from the data processing file accompanying:
###Lock, E. F., & Dunson, D. B. (2013). Bayesian consensus clustering. Bioinformatics (Oxford, England), 29(20), 2610–2616. http://doi.org/10.1093/bioinformatics/btt425

makeUMAPLabels <- function(umap_mat){
  
  x_med <- median(umap_mat[, 1])
  y_med <- median(umap_mat[, 2])
  
  my_labels <- rep(0, nrow(umap_mat))
  names(my_labels) <- row.names(umap_mat)
  
  my_labels[which(umap_mat[,1] > x_med & umap_mat[,2] > y_med)] <- 1
  my_labels[which(umap_mat[,1] > x_med & umap_mat[,2] <= y_med)] <- 2
  my_labels[which(umap_mat[,1] <= x_med & umap_mat[,2] <= y_med)] <- 3
  my_labels[which(umap_mat[,1] <= x_med & umap_mat[,2] > y_med)] <- 4

  my_labels <- as.factor(my_labels)
  
  my_labels
}

makeUMAPPlotData <- function(umap_mat, labels){
  
  plt_data <- data.frame(x = umap_mat[,1],
                         y = umap_mat[,2],
                         labels = labels
  )
  
  plt_data
}

plotUMAP <- function(plt_data){
  
  p <- ggplot(data = plt_data, aes(x = x, y = y, colour = labels)) +
    geom_point()
  
  p
}

# setwd("~/Dropbox/GwenDataIntegration/BCC")
# source("https://bioconductor.org/biocLite.R")

# if (!requireNamespace("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")

# BiocManager::install("impute")

library(impute)
library("MCMCpack")
library(magrittr)
library(stringr)
library(umap)
library(ggfortify)
library(ggplot2)

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


rowSums(miRNA.mat==0)
##Remove miRNAs with > 50% missing values


# === Analyse pre-processing ===================================================

# Create file names
save_loc <- "~/Documents/PhD/Year_1/Consensus_clustering/Analysis/BCC_TCGA_data/Data/"
file_type <- ".png"
plots <- c("_pre", "_post")

data_types <- c(
  "gene_expression",
  "methylation",
  "miRNA",
  "protein"
)

file_names <- lapply(plots, function(x){
  paste0(save_loc, data_types, x, file_type)
})

umap_file_names <- paste0(save_loc, data_types, "_umap", file_type)

# === Gene expression pre-processing ===========================================

# Get the Gene expression data in the right format and remove NAs
Exp.mat      <- as.matrix(GE[,2:349])
Exp.mat      <- impute.knn(Exp.mat) ##Impute missing values via KNN (K=10)
Exp.mat      <- Exp.mat$data %>% 
  set_rownames(GE[,1])

# Do PCA
ge_pca <- prcomp(Exp.mat)

# Visaulise the first two components
autoplot(ge_pca)

# Apply UMAP to the data
ge_umap <- umap(Exp.mat)

# Find labels from the UMAP to track how transform changes layout of data
ge_labels <- makeUMAPLabels(ge_umap$layout)
ge_plt_data <- makeUMAPPlotData(ge_umap$layout, ge_labels)
plotUMAP(ge_plt_data)  +
  labs(title = "Gene expression: pre-processing UMAP",
       subtitle = "Coloured by UMAP coordinates",
       x = "UMAP 1",
       y = "UMAP 2") 

# Save UMAP plot
ggsave(umap_file_names[1])

# Plot PCA with UMAP defined labels
autoplot(ge_pca, data = Exp.mat, colour = ge_labels) +
  labs(title = "Gene expression: pre-processing",
       subtitle = "Coloured by UMAP coordinates")

ggsave(file_names[[1]][1])

processedExpression  <- Exp.mat[apply(Exp.mat,1,sd)>1.5,] ###Filter to select only most variable genes

# Look at PCs after transform with the same labelling
p_ge_pca <- prcomp(processedExpression)
autoplot(p_ge_pca)
autoplot(p_ge_pca, data = processedExpression, colour = ge_labels[apply(Exp.mat,1,sd)>1.5]) +
  labs(title = "Gene expression: post-processing",
       subtitle = "Coloured by UMAP coordinates")

ggsave(file_names[[2]][1])

p_ge_umap <- umap(processedExpression)
p_ge_plt_data <- makeUMAPPlotData(p_ge_umap$layout, ge_labels[apply(Exp.mat,1,sd)>1.5])
plotUMAP(p_ge_plt_data)  +
  labs(title = "Gene expression: post-processing UMAP",
       subtitle = "Coloured by UMAP coordinates",
       x = "UMAP 1",
       y = "UMAP 2") 

ggsave("~/Documents/PhD/Year_1/Consensus_clustering/Analysis/BCC_TCGA_data/Data/gene_expression_umap_post_processing.png")

# === Methylation pre-processing ===============================================

meth_pca <- prcomp(Meth.mat)
autoplot(meth_pca)
meth_umap <- umap(Meth.mat)
plot(meth_umap$layout)

meth_labels <- makeUMAPLabels(meth_umap$layout)
meth_plt_data <- makeUMAPPlotData(meth_umap$layout, meth_labels)
plotUMAP(meth_plt_data)  +
  labs(title = "Methylation: pre-processing UMAP",
       subtitle = "Coloured by UMAP coordinates",
       x = "UMAP 1",
       y = "UMAP 2") 

ggsave(umap_file_names[2])


autoplot(meth_pca, data = Meth.mat, colour = meth_labels) +
  labs(title = "Methylation: pre-processing",
       subtitle = "Coloured by UMAP coordinates")

ggsave(file_names[[1]][2])


processedMethylation <- sqrt(Meth.mat)    ##take square root of methylation data
p_meth_pca <- prcomp(processedMethylation)

autoplot(p_meth_pca, data = processedMethylation, colour = meth_labels) +
  labs(title = "Methylation: post-processing",
       subtitle = "Coloured by UMAP coordinates")

ggsave(file_names[[2]][2])


p_meth_umap <- umap(processedMethylation)
p_meth_plt_data <- makeUMAPPlotData(p_meth_umap$layout, meth_labels)
plotUMAP(p_meth_plt_data) +
  labs(y="Petal length (cm)", x = "Sepal length (cm)")
  

# === miRNA pre-processing =====================================================

# Apply PCA for global structure
miRNA_pca <- prcomp(miRNA.mat)

# Apply UMAP to see local structure (should not be strongly affected by the transforms)
miRNA_umap <- umap(miRNA.mat)

# Create labels to record apporximate clustering from UMAP
miRNA_labels <- makeUMAPLabels(miRNA_umap$layout)

# Create a data.frame of UMAP coordinates and labelling
miRNA_plt_data <- makeUMAPPlotData(miRNA_umap$layout, miRNA_labels)

# Plot
plotUMAP(miRNA_plt_data) +
  labs(title = "miRNA: pre-processing UMAP",
       subtitle = "Coloured by UMAP coordinates",
       x = "UMAP 1",
       y = "UMAP 2") 

ggsave(umap_file_names[3])


# miRNA_exclude <- row.names(miRNA_plt_data)[miRNA_plt_data$x >10]
# autoplot(prcomp(miRNA.mat[! row.names(miRNA.mat) %in% miRNA_exclude, ]))
autoplot(miRNA_pca, data = miRNA.mat, colour = miRNA_labels, scale = 0) +
  labs(title = "miRNA: pre-processing",
       subtitle = "Coloured by UMAP coordinates")  # +
  # xlim(4e4, 5e4) +
  # ylim(5.5e3, 6.1e3)
# miRNA_pca$x[(miRNA_pca$x[,1] < 4e4),1:2]
# 
# summary(miRNA_pca$x[,1:2])
# 
# miRNA_exclude <- row.names(miRNA_pca$x)[(miRNA_pca$x[,1] < 4e4 & miRNA_pca$x[,2] < 6.1e3)]
# indices_to_drop <- ! row.names(miRNA.mat) %in% miRNA_exclude
# miRNA_reduced <- miRNA.mat[which(indices_to_drop),]
# miRNA_labels_reduced <- miRNA_labels[which(indices_to_drop)]
# 
# miRNA_pca_red <- prcomp(miRNA_reduced)
# 
# autoplot(miRNA_pca_red, data = miRNA_reduced, colour = miRNA_labels_reduced) +
#   labs(title = "miRNA: pre-processing",
#        subtitle = "Coloured by UMAP coordinates") 

# Save the PCA plot with UMAP colouring
ggsave(file_names[[1]][3])


# Transform the data
miRNA_to_drop <- rowSums(miRNA.mat==0) < 348*0.5
miRNA.mat    <- miRNA.mat[which(miRNA_to_drop),]
processedmiRNA       <- log(1+miRNA.mat) ##take log of miRNA data

# Take the PCA of the transformed data
p_miRNA_pca <- prcomp(processedmiRNA)

# Plot with the same labelling as previosuly applied
autoplot(p_miRNA_pca, data = processedmiRNA, colour = miRNA_labels[which(miRNA_to_drop)]) +
  labs(title = "miRNA: post-processing",
       subtitle = "Coloured by UMAP coordinates") 

ggsave(file_names[[2]][3])

p_miRNA_umap <- umap(processedmiRNA)
p_miRNA_plt_data <- makeUMAPPlotData(p_miRNA_umap$layout, miRNA_labels[which(miRNA_to_drop)])
plotUMAP(p_miRNA_plt_data) +
  labs(title = "miRNA: post-processing UMAP",
       subtitle = "Coloured by UMAP coordinates",
       x = "UMAP 1",
       y = "UMAP 2") 

ggsave("~/Documents/PhD/Year_1/Consensus_clustering/Analysis/BCC_TCGA_data/Data/miRNA_umap_post_processing.png")


# === Protein pre-processing ===================================================

protein_pca <- prcomp(Protein.mat)
autoplot(protein_pca)
protein_umap <- umap(Protein.mat)
plot(protein_umap$layout)

protein_labels <- makeUMAPLabels(protein_umap$layout)
protein_plt_data <- makeUMAPPlotData(protein_umap$layout, protein_labels)
plotUMAP(protein_plt_data) +
  labs(title = "Protein: pre-processing UMAP",
       subtitle = "Coloured by UMAP coordinates",
       x = "UMAP 1",
       y = "UMAP 2") 

ggsave(umap_file_names[4])

autoplot(protein_pca, data = Protein.mat, colour = protein_labels) +
  labs(title = "Protein: pre-processing",
       subtitle = "Coloured by UMAP coordinates") 

ggsave(file_names[[1]][4])


processedProtein     <- scale(Protein.mat,center=TRUE,scale=TRUE) #Column center/scale protein
p_protein_pca <- prcomp(processedProtein)
autoplot(p_protein_pca, data = processedProtein, colour = protein_labels)  +
  labs(title = "Protein: post-processing",
       subtitle = "Coloured by UMAP coordinates") 

ggsave(file_names[[2]][4])


p_protein_umap <- umap(processedProtein)
p_protein_plt_data <- makeUMAPPlotData(p_protein_umap$layout, protein_labels)
plotUMAP(p_protein_plt_data) +
  labs(title = "Protein: post-processing UMAP",
       subtitle = "Coloured by UMAP coordinates",
       x = "UMAP 1",
       y = "UMAP 2") 

# === MDI pre-processing =======================================================

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
