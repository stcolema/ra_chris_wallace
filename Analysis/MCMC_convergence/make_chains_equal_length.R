#!/usr/env/bin/Rscript

# Find the shortest chain and cut any iterations beyond this for all chains 
# (important for tests of convergence across chains)

library(optparse)

input_arguments <- function() {
  option_list <- list(
    
    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("-d", "--dir"),
                          type = "character",
                          default = "./",
                          help = "Directory containing files to be reduced [default= %default]",
                          metavar = "character"
    ),
    
    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("-n", "--n_iter"),
                          type = "integer",
                          default = NA,
                          help = "Number of iterations to keep in each chain [default= %default]",
                          metavar = "integer"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

args <- input_arguments()
data_dir <- args$dir
smallest <- args$n_iter

if(is.na(smallest)){
  smallest <- 1e8
}

# The main directory
# data_dir <- "~/Documents/PhD/Year_1/Consensus_clustering/Analysis/MDI_test_data/Easier_convergence/MATLAB_MDI/N_clust_50/"
dir.create(paste0(data_dir, "Original_data"))
dir.create(paste0(data_dir, "Abbreviated_data"))

files <- list.files(data_dir, include.dirs = F, pattern = ".csv")
n_files <- length(files)
my_files <- list()
# smallest <- 1e7
for(i in 1:n_files){
  my_files[[i]] <- .curr_data <- read.csv(paste0(data_dir, files[i]))
  
  smallest <- min(nrow(.curr_data), smallest)
  
  
}


my_cut_files <- list()
for(i in 1:n_files){
  my_cut_files <- .small_data <- my_files[[i]][1:smallest, ]
  write.csv(.small_data, paste0(data_dir, "Abbreviated_data/", files[i]), row.names = F)
}

