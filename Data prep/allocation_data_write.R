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

# === Libraries ================================================================

source("~/Desktop/MDI/mdipp-1.0.1/scripts/analysis.R")
library(magrittr)
library(dplyr)
library(data.table)
library(rlist)
library(pheatmap) # install.packages("pheatmap", dep = T)
library(RColorBrewer)
# output <- loadDataGauss("~/Desktop/First attempt/output_1.csv")

# === Setup ====================================================================

setwd("~/Desktop/MDI/Data/Full_sets")

# Read in the probes present - this is a matrix of bools with the first column
# as the probe IDs and the remaining columns corrresponding to the cell types
# with TRUE indicating the probe is present in this cell type (i.e. not added
# manually with an imputed value) and FALSE indicates we added it in.
probes_present_dt <- fread("~/Desktop/MDI/Data/probes_present_per_dataset.csv")

# Read in the file relating the probe IDs to  the related gene
probe_key <- fread("~/Desktop/MDI/Data/probe_key.csv")

# Read in the MDI output file
mdi_output_file <- "output_full_sets.csv"

# Convert this into useable output using functions provided by Sam Mason
mcmc_output <- readMdiMcmcOutput(mdi_output_file)

# Declare empty variable to capture information
n_genes <- NA
mdi_allocation <- list()
allocation_list <- list()

# MDI call specific values
num_datasets <- 9
n_iter <- 1e4
thin <- 25
burn <- 0.1 * (n_iter / thin)
eff_n_iter <- n_iter / thin - burn

# Create a dataframe to hold the output
col_names <- paste0("D", 1:num_datasets)

# The number of genes is the number of columns in the mcmc_output
# exclusing the columns containing information on the phi and mass parameters
if (is.na(n_genes)) {
  n_genes <- mcmc_output %>%
    select(contains("Dataset")) %>%
    ncol() / num_datasets
}

# Create an empty dataframe with column names corresponding to dataset numbers
compare_df <- data.frame(matrix(ncol = num_datasets, nrow = n_genes))
# colnames(compare_df) <- col_names

# The actual names of the datasets are
files_present <- list.files(path = ".") %>%
  grep(".csv", ., value = TRUE)

# Remove the output file from this list if you saved it in the same location as
# the data files
if (mdi_output_file %in% files_present) {
  index_to_empty <- files_present == mdi_output_file
  files_present %<>%
    list.remove(index_to_empty)
}
# Remove the file extension
dataset_names <- tools::file_path_sans_ext(files_present)

# Have to have the order from the call
# CD4, CD8, CVD19, CD14, CD15, platelets, ileonic, colonic and rectal biopsies
# transcriptome data for six circulating immune cell types (CD4+ T lymphocytes,
# CD8+ T lymphocytes, CD19+ B lymphocytes, CD14+ monocytes, CD15+ granulocytes,
# platelets) as well as ileal, colonic, and rectal biopsies (IL, TR, RE)
# 323 healthy Europeans
# files_present <- c("CD14.csv", "CD15.csv", "IL.csv", "CD4.csv", "CD8.csv", "CD19.csv", "PLA.csv", "RE.csv", "TR.csv")
dataset_names <- tools::file_path_sans_ext(files_present)

# Set the cell type to the column names (make sure order is as per MDI call)
colnames(compare_df) <- dataset_names

# === MDI output ===============================================================

# Capture the allocation information in the named lists and the predicted
# allocation in the dataframe
for (i in 1:num_datasets) {
  dataset_name <- paste0("D", i)
  
  # Get the allocation and drop the burn in
  mdi_allocation[[i]] <- .mdi_alloc <- getMdiAllocations(mcmc_output, i) %>%
    .[-(1:burn), ]
  
  allocation_list[[i]] <- .pred_alloc <- apply(.mdi_alloc, 2, median)
  
  compare_df[[dataset_names[i]]] <- .pred_alloc
  
  if (i == 1) {
    row.names(compare_df) <- names(.pred_alloc)
  }
}

# Possibly need to rearrange based on call order
# compare_df <- compare_df[, c(4,5,1,2,6,7,3,8,9)]

# Check the number and vlaues of clusters present
clusters_present <- apply(compare_df, 2, unique)

# Set the row names
row.names(probes_present_dt) <- probes_present_dt$V1

# Remove the ID column
probes_present_dt_no_V1 <- probes_present_dt[, -1]

# Create a new dataframe with 0's for the probes not present in that dataframe
# and the allocation label otherwise
compare_df_2 <- as.matrix(compare_df) * as.matrix(probes_present_dt_no_V1)

# compare_df_2[compare_df_2 == 0] <- NA

# Add the row and column names
row.names(compare_df_2) <- row.names(compare_df)
colnames(compare_df_2) <- colnames(compare_df)

# Check if this looks correct
head(compare_df_2)

# Move to a new labelling system (as MDI overfits the number of clusters we have
# an unevenly distributed range of values - in our new system there is equal
# distance between each label)
# sort(unique(c(as.matrix(compare_df))))

n_clust <- length(unique(c(compare_df_2)))
old_labels <- sort(unique(c(compare_df_2)))
new_labels <- 1:n_clust
new_labels[1] <- -1
key <- data.frame(old = old_labels, new = new_labels)

# Create a dataframe with the new labels
compare_df_new_labels <- compare_df_2

# Use [] to ensure that the structure is preserved
compare_df_new_labels[] <- key$new[match(unlist(compare_df_2), key$old)]


# === Heatmapping ==============================================================

col_pal <- col_pal_old

# Pull out the Gene IDs in the correct order
gene_id <- probe_key %>%
  .[match(.$ProbeID, row.names(compare_df_new_labels))] %>%
  .$Gene

# Move from Probe ID to Gene ID
row.names(compare_df_new_labels) <- gene_id

compare_df_new_labels <- as.data.table(compare_df_new_labels)

compare_df_new_labels$V1 <- gene_id

fwrite(compare_df_new_labels_sc, "allocation _data_run_2.csv")
