#!/usr/env/bin Rscript

# /home/MINTS/sdc56/Desktop/MDI/mdipp-1.0.1/mdipp
# N CD14.csv
# N CD15.csv
# N CD19.csv
# N CD4.csv
# N CD8.csv
# N IL.csv
# N PLA.csv
# RE.csv
# TR.csv
# -n 10000
# -t 25
# -s 1
# > output_full_sets.csv


# Function to find the mode of a vector
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

setwd("/home/MINTS/sdc56/Desktop/ra_chris_wallace")

# === Libraries ================================================================

source("~/Desktop/MDI/mdipp-1.0.1/scripts/analysis.R") # install.packages("mcclust", dep = T)
Rcpp::sourceCpp("Analysis/posterior_sim_mat.cpp") # install.packages("Rcpp", dep = T)
library(magrittr)
library(tibble) # for dataframe of lists
library(dplyr)
library(data.table)
# library(rlist) # install.packages("rlist", dep = T)
library(pheatmap) # install.packages("pheatmap", dep = T)
library(RColorBrewer)
library(mclust) # install.packages("mclust", dep = T)
library(ggplot2)

# output <- loadDataGauss("~/Desktop/First attempt/output_1.csv")

# === Setup ====================================================================

# Have we already done the initial mdi output analysis and prep
data_written <- F

# Instructions to check, save and do several different things
save_plots <- T
plot_type <- ".pdf"
do_heatplot_clusterings <- T
do_rand_plot <- T
do_individual_comparison <- T

# setwd("~/Desktop/MDI/Data/Full_sets")

# Read in the probes present - this is a matrix of bools with the first column
# as the probe IDs and the remaining columns corrresponding to the cell types
# with TRUE indicating the probe is present in this cell type (i.e. not added
# manually with an imputed value) and FALSE indicates we added it in.
probes_present_dt <- fread("Analysis/probes_present_per_dataset.csv")

# Read in the file relating the probe IDs to  the related gene
probe_key <- fread("Analysis/probe_key.csv")

# Read in the MDI output file
file_path <- "Analysis/MDI_runs/vsn_many_seeds/"
mdi_output_files <- list.files(path = file_path, full.names = T, include.dirs = F) %>% 
  grep("csv", ., value = TRUE)

num_files <- length(mdi_output_files)
file_names <- basename(tools::file_path_sans_ext(mdi_output_files))

# Declare empty variable to capture information
n_genes <- NA
mdi_allocation <- list()
allocation_list <- list()

# MDI call specific values
num_datasets <- 9
n_iter <- 800
thin <- 1
burn <- 0

eff_n_iter <- n_iter / thin # - burn

# Create a dataframe to hold the output
col_names <- paste0("D", 1:num_datasets)

# The number of genes is the number of columns in the mcmc_output
# exclusing the columns containing information on the phi and mass parameters
n_genes <- 18524 # SPECIFIC TO CURRENT RUN


# The actual names of the datasets are
# files_present <- list.files(path = "~/Desktop/MDI/Data/VSN_NA_data") %>%
# grep(".csv", ., value = TRUE)
files_present <- c(
  "CD14.csv",
  "CD15.csv",
  "CD19.csv",
  "CD4.csv",
  "CD8.csv",
  "IL.csv",
  "PLA.csv",
  "RE.csv",
  "TR.csv"
)

# Remove the file extension
dataset_names <- tools::file_path_sans_ext(files_present)

# mdi_output_file <- "Analysis/MDI_runs/vsn_many_seeds/out_seed_1.csv"
# If already found allocations and written to file, re-load
if (data_written) {
  compare_tibble <- readRDS("~/Desktop/ra_chris_wallace/Analysis/MDI_runs/vsn_many_seeds/compare_tibble.rds")
}


if (!data_written) {
  # Convert this into useable output using functions provided by Sam Mason
  mcmc_out_lst <- list() # vector("list", num_files)
  
  for (i in 1:num_files) {
    curr_name <- file_names[[i]]
    mcmc_out_lst[[curr_name]] <- readMdiMcmcOutput(mdi_output_files[[i]])
  }
  
  
  # mcmc_output <- readMdiMcmcOutput(mdi_output_file)
  
  clust_occ <- mcmc_out_lst %>%
    lapply(getClustersOccupied)
  
  # clust_occ <- getClustersOccupied(mcmc_output)
  # head(clust_occ)
  
  # By the 8 x 25 = 200 iteration we have reached stationarity in the number of
  # clusters occupied
  # summary(clust_occ[8:800, ])
  
  # NOte that as each mdi output is on the same data, we only need to count one dataset
  if (is.na(n_genes)) {
    n_genes <- mcmc_out_lst[[1]] %>%
      select(contains("Dataset")) %>%
      ncol() / num_datasets
  }
  
  
  
  # Create an empty dataframe with column names corresponding to dataset numbers
  compare_df <- data.frame(matrix(ncol = num_datasets, nrow = n_genes))
  # colnames(compare_df) <- col_names
  
  # Remove the output file from this list if you saved it in the same location as
  # the data files
  if (mdi_output_file %in% files_present) {
    index_to_empty <- files_present == mdi_output_file
    files_present %<>%
      list.remove(index_to_empty)
  }
  
  # Have to have the order from the call
  # CD4, CD8, CVD19, CD14, CD15, platelets, ileonic, colonic and rectal biopsies
  # transcriptome data for six circulating immune cell types (CD4+ T lymphocytes,
  # CD8+ T lymphocytes, CD19+ B lymphocytes, CD14+ monocytes, CD15+ granulocytes,
  # platelets) as well as ileal, colonic, and rectal biopsies (IL, TR, RE)
  # 323 healthy Europeans
  # files_present <- c("CD14.csv", "CD15.csv", "IL.csv", "CD4.csv", "CD8.csv", "CD19.csv", "PLA.csv", "RE.csv", "TR.csv")
  # dataset_names <- tools::file_path_sans_ext(files_present)
  
  # Set the cell type to the column names (make sure order is as per MDI call)
  colnames(compare_df) <- dataset_names
  
  # Now put everything in a tibble
  compare_tibble <- tibble(
    mdi = rep(file_names, num_datasets),
    dataset = rep(dataset_names, num_files), # unlist(lapply(dataset_names, rep, num_files))
    seed = unlist(lapply(1:num_files, rep, num_datasets)), # rep(1:num_files, num_datasets), 
    n_genes = n_genes,
    pred_allocation = list(compare_df),
    mdi_allocation = rep(vector("list", num_files), num_datasets)
  )
  
  # === MDI output ===============================================================
  
  # x4 <- sample(1:eff_iter, 50, replace=F)
  
  for (j in 1:num_files) {
    # compare_tibble$MDI_allocation[j] <- list()
    # compare_tibble$Pred_allocation[j] <- list()
    
    mdi_pos_sim_mat <- list()
    check_median_makes_sense_map <- list()
    
    # Capture the allocation information in the named lists and the predicted
    # allocation in the dataframe
    
    
    for (i in 1:num_datasets) {
      dataset_name <- paste0("D", i)
      
      # Get the allocation and drop the burn in
      mdi_allocation[[i]] <- .mdi_alloc <- getMdiAllocations(mcmc_out_lst[[j]], i)
      compare_tibble$mdi_allocation[i + (j - 1) * num_files] <- .mdi_alloc
      
      # sams_pos_sim <- genPosteriourSimilarityMatrix(.mdi_alloc)
      #
      # # Individuals are columns rather than rows (henc transpose)
      # mdi_pos_sim_mat[[i]] <- list()
      # mdi_pos_sim_mat[[i]][[1]] <- .sim_mat <- similarity_mat(t(.mdi_alloc[,1:1000]))
      # mdi_pos_sim_mat[[i]][[2]] <- .sim_mat <- similarity_mat(t(.mdi_alloc[,5000:6000]))
      # mdi_pos_sim_mat[[i]][[3]] <- .sim_mat <- similarity_mat(t(.mdi_alloc[,11000:12000 ]))
      # # break
      # # By checking the imilarity of each row we can decide if the median is an
      # # accurate method to allocate class (if all rows are highly similar label
      # # flipping did not occur)
      # check_median_makes_sense_map[[i]] <- .sense_check <- similarity_mat(.mdi_alloc[c(seq(1, eff_n_iter, 25), eff_n_iter), ])
      # break
      # }
      
      # pheatmap(mdi_pos_sim_mat[[1]])
      
      allocation_list[[i]] <- .pred_alloc <- apply(.mdi_alloc[-(1:burn), ], 2, getmode)
      # compare_df[[dataset_names[i]]] <- .pred_alloc
      
      compare_tibble$pred_allocation[i + (j - 1) * num_files] <- .pred_alloc
      
      
      # if (i == 1) {
      # row.names(compare_df) <- names(.pred_alloc)
      # }
    }
    # compare_tibble$MDI_allocation[j] <- mdi_allocation
    # compare_tibble$Pred_allocation[j] <- compare_df
  }
  
  saveRDS(
    compare_tibble,
    "~/Desktop/ra_chris_wallace/Analysis/MDI_runs/vsn_many_seeds/compare_tibble.rds"
  )
}