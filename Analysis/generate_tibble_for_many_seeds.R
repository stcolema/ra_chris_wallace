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



# === Libraries ================================================================

source("~/Desktop/MDI/mdipp-1.0.1/scripts/analysis.R") # install.packages("mcclust", dep = T)
Rcpp::sourceCpp("~/Desktop/ra_chris_wallace/Analysis/posterior_sim_mat.cpp") # install.packages("Rcpp", dep = T)
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
save_plots <- F
plot_type <- ".png"
do_heatplot_clusterings <- T
do_rand_plot <- T
do_individual_comparison <- T

# setwd("~/Desktop/MDI/Data/Full_sets")

# Read in the probes present - this is a matrix of bools with the first column
# as the probe IDs and the remaining columns corrresponding to the cell types
# with TRUE indicating the probe is present in this cell type (i.e. not added
# manually with an imputed value) and FALSE indicates we added it in.

# Generic base for the data
gen_base <- "/home/MINTS/sdc56/Desktop/Many seeds/"

# The gene sets
gene_sets <- c("All_seeds/")
n_sets <- length(gene_sets)

# Create the path to each of the gene sets data
gene_set_wd <- gen_base %>%
  paste0(gene_sets)

# Meta data directory name in each of the gene set directories
meta_data_wd <- "Meta_data"

# The name of the files containing the information about which probes were
# actually present in each dataset and which were added as a row of 0's
probes_present_csv <- "probes_present_per_dataset.csv"

# Empty list to hold information about each gene set
probes_present_dt <- vector("list", n_sets)
names(probes_present_dt) <- gene_sets

# Generic file path within gene set directories for MDI output
mdi_out_file_path <- "hpc_output"

# Read in the MDI output file
file_path <- gene_set_wd %>%
  paste0(mdi_out_file_path)

# List to hold names of MDI output files
mdi_output_files <- vector("list", n_sets)
names(mdi_output_files) <- gene_sets

# Record the number of files in each gene set
num_files <- vector("list", n_sets)
names(num_files) <- gene_sets

# List to record the file names without extensions
file_names <- vector("list", n_sets) #
names(file_names) <- gene_sets

# Iterate through the sets and read in the information about which probes are in
# which datasets
for (i in 1:n_sets) {
  # Record the current set
  curr_set <- gene_sets[[i]]

  # Create the file name to read in
  .f_i <- paste(gene_set_wd[[i]], meta_data_wd, probes_present_csv, sep = "/")

  # Read in the file relating the probe IDs to  the related gene
  probes_present_dt[[curr_set]] <- fread(.f_i)

  # Read in the files present in the directories
  mdi_output_files[[curr_set]] <- list.files(path = file_path[[i]], full.names = T, include.dirs = F) %>%
    grep("csv", ., value = TRUE)

  # The stripped file names
  file_names[[curr_set]] <- mdi_output_files[[curr_set]] %>%
    tools::file_path_sans_ext() %>%
    basename()

  # This is for handling the fact that 10 comes before 9 in strings but not in
  # numbers. Possibly a bad idea.
  
  seeds_present <- file_names[[curr_set]] %>%
    stringr::str_extract(., "\\-*\\d+\\.*\\d*") %>%
    as.numeric()
  
  file_name_order <- seeds_present %>%
    order()

  # Re-order
  mdi_output_files[[curr_set]] <- mdi_output_files[[curr_set]][file_name_order]
  file_names[[curr_set]] <- file_names[[curr_set]][file_name_order]

  # The number of otuput files for each gene set MDI
  num_files[[curr_set]] <- mdi_output_files[[curr_set]] %>%
    length()
}

# Declare empty variable to capture information
n_genes <- NA
mdi_allocation <- list()
allocation_list <- list()

# MDI call specific values
num_datasets <- 9
n_iter <- 500
thin <- 500
burn <- 0

eff_n_iter <- n_iter / thin - burn

# Create a dataframe to hold the output
col_names <- paste0("D", 1:num_datasets)

# The number of genes is the number of columns in the mcmc_output
# exclusing the columns containing information on the phi and mass parameters
# n_genes <- 18524 # SPECIFIC TO CURRENT RUN


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

# === Massive for loop =========================================================

save_name <- paste(file_path, "compare_tibble.rds", sep = "/")

gene_sets_out <- vector("list", n_sets)
names(gene_sets_out) <- gene_sets

# mdi_output_file <- "Analysis/MDI_runs/vsn_many_seeds/out_seed_1.csv"
# If already found allocations and written to file, re-load
if (data_written) {
  for (i in 1:n_sets) {
    curr_file <- save_name[[i]]
    set_name <- gene_sets[[i]]
    curr_tibble <- readRDS(curr_file)
    gene_sets_out[[set_name]] <- curr_tibble
  }
  # compare_tibble <- readRDS("~/Desktop/ra_chris_wallace/Analysis/MDI_runs/vsn_many_seeds/compare_tibble.rds")
}


if (!data_written) {
  n_genes <- vector("list", n_sets)
  names(n_genes) <- gene_sets

  for (l in 1:n_sets) {
    curr_set <- gene_sets[[l]]

    # Convert this into useable output using functions provided by Sam Mason
    mcmc_out_lst <- list() # vector("list", num_files)

    # Read in the MDI output for each file
    for (i in 1:num_files[[curr_set]]) {
    # for (i in 1:2) {
      curr_name <- file_names[[curr_set]][[i]]
      mcmc_out_lst[[curr_name]] <- readMdiMcmcOutput(mdi_output_files[[curr_set]][[i]])
    }


    # mcmc_output <- readMdiMcmcOutput(mdi_output_file)

    # Record the number of occpied clusters
    clust_occ <- mcmc_out_lst[[curr_set]] %>%
      lapply(getClustersOccupied)


    # NOte that as each mdi output is on the same data, we only need to count one dataset
    if (is.null(n_genes[[curr_set]])) {
      n_genes[[curr_set]] <- .n_gene_curr <- mcmc_out_lst[[1]] %>%
        select(contains("Dataset")) %>%
        ncol() / num_datasets
    }

    # Matrix to hold allocation values
    alloc_matrix <- data.frame(matrix(nrow = eff_n_iter, ncol = n_genes[[curr_set]]))

    # Create an empty dataframe with column names corresponding to dataset numbers
    compare_df <- data.frame(matrix(ncol = num_datasets, nrow = n_genes[[curr_set]]))
    # colnames(compare_df) <- col_names

    # Remove the output file from this list if you saved it in the same location as
    # the data files
    # if (mdi_output_file %in% files_present) {
    #   index_to_empty <- files_present == mdi_output_file
    #   files_present %<>%
    #     list.remove(index_to_empty)
    # }

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

    # Column containing the random seed for the tibble
    seed_var <- 1:num_files[[curr_set]] %>%
      lapply(., rep, num_datasets) %>%
      unlist()

    # Column contatining the tissue type
    tissue_var <- dataset_names %>%
      rep(num_files[[curr_set]])

    # The mdi output file name
    f_name <- file_names[[curr_set]] %>% lapply(rep, num_datasets) %>% unlist() %>% sort()
    assoc_seeds <- f_name %>% stringr::str_extract(., "\\-*\\d+\\.*\\d*") %>% as.numeric()

    f_name <- f_name[order(assoc_seeds)]

    # Now put everything in a tibble
    compare_tibble <- tibble(
      mdi = f_name,
      dataset = tissue_var, # unlist(lapply(dataset_names, rep, num_files))
      seed = seed_var, # rep(1:num_files, num_datasets),
      n_genes = .n_gene_curr,
      # pred_allocation = rep(vector("list", num_files[[curr_set]]), num_datasets), # list(compare_df),
      mdi_allocation = list(alloc_matrix) # rep(vector("list", num_files[[curr_set]]), num_datasets)
    )

    # === MDI output ===============================================================

    # x4 <- sample(1:eff_iter, 50, replace=F)
    # for (j in 1:2) {
    for (j in 1:num_files[[curr_set]]) {
      # compare_tibble$MDI_allocation[j] <- list()
      # compare_tibble$Pred_allocation[j] <- list()

      curr_seed <- seeds_present[[j]]
      
      mdi_pos_sim_mat <- list()
      check_median_makes_sense_map <- list()

      # Capture the allocation information in the named lists and the predicted
      # allocation in the dataframe


      for (i in 1:num_datasets) {
        dataset_name <- dataset_names[[i]]

        tibble_row <- row.names(compare_tibble) %>%
          .[compare_tibble$dataset == dataset_name & compare_tibble$seed == curr_seed] %>%
          as.numeric()

        # Get the allocation and drop the burn in
        mdi_allocation[[i]] <- .mdi_alloc <- getMdiAllocations(mcmc_out_lst[[j]], i)
        compare_tibble$mdi_allocation[tibble_row][[1]] <- .mdi_alloc

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

        # if (burn > 0) {
        #   allocation_list[[i]] <- .pred_alloc <- apply(.mdi_alloc[-(1:burn), ], 2, getmode)
        # } else {
        #   allocation_list[[i]] <- .pred_alloc <- apply(.mdi_alloc, 2, getmode)
        # }
        # # compare_df[[dataset_names[i]]] <- .pred_alloc
        # 
        # compare_tibble$pred_allocation[tibble_row][[1]] <- .pred_alloc


        # if (i == 1) {
        # row.names(compare_df) <- names(.pred_alloc)
        # }
      }
      # compare_tibble$MDI_allocation[j] <- mdi_allocation
      # compare_tibble$Pred_allocation[j] <- compare_df
    }



    saveRDS(
      compare_tibble,
      save_name[[l]]
    )

    gene_sets_out[[curr_set]] <- compare_tibble
  }
}
