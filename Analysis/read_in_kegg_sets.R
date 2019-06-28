#!/usr/bin/env Rscript

# Function to read in KEGG sets

library(data.table)
library(magrittr)


find_gene_sets <- function(file_name, universe, unique_universe = NULL){
  
  # Read in the file
  f <- scan(file_name,"character()",sep="\n") %>% strsplit(.,"\t")
  
  # Create a list of gene sets with names based on pathways
  hm <- lapply(f,function(x) {
    
    # Remove the pathway name and the url
    genes <- x[c(-1,-2)] 
    
    # Ensure the genes are present in our universe of genes
    genes <- genes[genes %in% universe]
    
    # If given a unique identifier (pretty much for CEDAR dataset some probes 
    # map to the same gene so we have duplicates and thus require unique gene
    # identifiers)
    if(! is.null(unique_universe)){
      genes <- my_universe[universe %in% genes]
    }
    
    genes
  }) %>% set_names(sapply(f,'[[',1)) # set names as the pathway names
  
  hm
}

# probe_key <- fread("./Analysis/probe_key_unique_names.csv")
# universe <- probe_key$Gene
# my_universe <- probe_key$Unique_gene_name
# 
# 
# f <- scan("../kegg_msigdb.txt","character()",sep="\n") %>% strsplit(.,"\t")
# hm <- lapply(f,function(x) {
#   genes <- x[c(-1,-2)] 
#   genes <- genes[genes %in% universe]
#   re_genes <- my_universe[universe %in% genes]
# }) %>% set_names(sapply(f,'[[',1))