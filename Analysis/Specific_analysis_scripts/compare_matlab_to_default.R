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

# Function to plot and save tree for a given similarity matrix
plot_tree <- function(sim_mat, plot_name, plot_type = ".png"){
  
  if(! (plot_type %in% c(".pdf", ".png"))){
    stop("Wrong plot type attempted: must be one of '.pdf' or '.png'")
  } 
  
  # Open a file
  if(plot_type == ".png"){
    png(plot_name) 
  } else {
    pdf(plot_name)
  }
  # 2. Create a plot
  sim_mat %>%
    dist() %>% 
    hclust() %>%
    plot()
  
  # Close the pdf file
  dev.off() 
}

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


# Read in the probes present - this is a matrix of bools with the first column
# as the probe IDs and the remaining columns corrresponding to the cell types
# with TRUE indicating the probe is present in this cell type (i.e. not added
# manually with an imputed value) and FALSE indicates we added it in.
probes_present_dt <- fread("/home/MINTS/sdc56/Desktop/MDI_small_geneset_outputs/Meta_data/probes_present_per_dataset.csv")

# Read in the file relating the probe IDs to  the related gene
probe_key <- fread("Analysis/probe_key.csv")

# Read in the MDI output file
matlab_file_path <- "/home/MINTS/sdc56/Desktop/MDI_small_geneset_outputs/matlab_output/"
mdi_file_path <- "/home/MINTS/sdc56/Desktop/MDI_small_geneset_outputs/Many_seeds_small_7_output/"
# file_path <- "Analysis/MDI_runs/vsn_many_seeds/"

matlab_compare_tibble <- readRDS(paste0(matlab_file_path, "compare_tibble.rds"))
mdi_compare_tibble <- readRDS(paste0(mdi_file_path, "compare_tibble.rds"))

probe_names <- colnames(mdi_compare_tibble$mdi_allocation[1][[1]])

# The actual names of the datasets are
# files_present <- list.files(path = "~/Desktop/MDI/Data/VSN_NA_data") %>%
# grep(".csv", ., value = TRUE)
files_present <- c(
  "CD14.csv",
  "CD19.csv",
  "CD4.csv",
  "CD8.csv",
  "IL.csv",
  "RE.csv",
  "TR.csv"
)
num_datasets <- length(files_present)

# Create a dataframe to hold the output
col_names <- paste0("D", 1:num_datasets)

# The number of genes is the number of columns in the mcmc_output
# exclusing the columns containing information on the phi and mass parameters
n_genes <- 102 # NA # SPECIFIC TO CURRENT RUN

# Remove the file extension
dataset_names <- tools::file_path_sans_ext(files_present)

# mdi_output_file <- "Analysis/MDI_runs/vsn_many_seeds/out_seed_1.csv"
# If already found allocations and written to file, re-load

# === Probes present ===========================================================

# Find which probes are relevant from the full set
probes_actually_present_ind <- probes_present_dt %>% 
  magrittr::use_series("V1") %>% 
  magrittr::is_in(probe_names)

# Select these
probes_actually_present <- probes_present_dt[probes_actually_present_ind, ]

# Find the appropriate order
probes_order <- match(probe_names, probes_actually_present$V1)

# Order the probes so comparable to allocation data frame
probes_present_ordered <- probes_actually_present[probes_order, ] %>% 
  as.data.frame()

# Remove the irrelevant columns and set row names
probes_present_final <- probes_present_ordered %>% 
  magrittr::set_rownames(probes_actually_present$V1) %>% 
  magrittr::extract( , -c(1, 3, 8)) %>% 
  set_colnames(dataset_names)

apply(probes_present_final, 2, sum)

#  === Heatmapping ==============================================================

# Possibly need to rearrange based on call order
compare_df <- data.frame(matrix(ncol = num_datasets, nrow = n_genes))

colnames(compare_df) <- dataset_names
row.names(compare_df) <- probe_names

matlab_compare_df <- mdi_compare_df <- compare_df

for (set in dataset_names){
  matlab_compare_df[[set]] <- matlab_compare_tibble$pred_allocation[matlab_compare_tibble$dataset == set][[1]]
  mdi_compare_df[[set]] <- mdi_compare_tibble$pred_allocation[mdi_compare_tibble$dataset == set][[1]]
}

# Check the number and vlaues of clusters present
# clusters_present <- lapply(compare_tibble$pred_allocation, unique)
matlab_clusters_present <- apply(matlab_compare_df, 2, unique)
mdi_clusters_present <- apply(mdi_compare_df, 2, unique)

# Create a new dataframe with 0's for the probes not present in that dataframe
# and the allocation label otherwise
matlab_compare_df_2 <- as.matrix(matlab_compare_df) * as.matrix(probes_present_final) %>% 
  as.data.frame() %>% 
  magrittr::set_rownames(., probe_names)

# compare_df_2[compare_df_2 == 0] <- NA

# Add the row and column names
row.names(matlab_compare_df_2) <- row.names(compare_df)
colnames(matlab_compare_df_2) <- colnames(compare_df)

# Check if this looks correct
head(matlab_compare_df_2)

mdi_compare_df_2 <- as.matrix(mdi_compare_df) * as.matrix(probes_present_final) %>% 
  as.data.frame() %>% 
  magrittr::set_rownames(., probe_names)

# mdi_compare_df_2[mdi_compare_df_2 == 0] <- NA

# Add the row and column names
row.names(mdi_compare_df_2) <- row.names(compare_df)
colnames(mdi_compare_df_2) <- colnames(compare_df)

# Check if this looks correct
head(mdi_compare_df_2)


# Move to a new labelling system (as MDI overfits the number of clusters we have
# an unevenly distributed range of values - in our new system there is equal
# distance between each label)
# sort(unique(c(as.matrix(compare_df))))

matlab_n_clust <- length(unique(c(matlab_compare_df_2)))
mdi_n_clust <- length(unique(c(mdi_compare_df_2)))

n_clust <- max(matlab_n_clust, mdi_n_clust)

col_pal_old <- c("#DDDDDD", rev(brewer.pal(n = n_clust, "RdYlBu"))) # name = "RdYlBu")))
col_pal <- c(rev(brewer.pal(n = n_clust, "RdYlBu"))) # name = "RdYlBu")))


col_pal <- col_pal_old

probe_key_rel <- probe_key[probe_key$ProbeID %in% probe_names,]

# Pull out the Gene IDs in the correct order
gene_id <- probe_key_rel %>%
  .[match(probe_names, .$ProbeID)] %>%
  .$Gene 

# Move from Probe ID to Gene ID
matlab_compare_df_genes <- as.matrix(matlab_compare_df_2) %>% 
  magrittr::set_rownames(gene_id)
mdi_compare_df_genes <- as.matrix(mdi_compare_df_2) %>% 
  magrittr::set_rownames(gene_id)


# Heatmap of allocaiton

# We check if the NA rows (filled with 0's) are all in the same cluster and alone in this cluster
ph_matlab_no_zeros <- pheatmap(matlab_compare_df,
                      cluster_cols = T,
                      color = col_pal,
                      cellheight = 3,
                      cellwidth = 60,
                      main = "MATLAB output with normalisation"
)

zerod_row_order <- ph_matlab_no_zeros$tree_row[["order"]]
zerod_col_order <- ph_matlab_no_zeros$tree_col[["order"]]

matlab_zero_comparison <- matlab_compare_df_genes[zerod_row_order, zerod_col_order]

ph_matlab_zerod <- pheatmap(matlab_zero_comparison,
                      cluster_cols = F,
                      cluster_rows = F,
                      color = col_pal,
                      cellheight = 3,
                      cellwidth = 60,
                      main = "MATLAB output with normalisation - NAs in 0th cluster"
)
# They are very much not

ph_matlab <- pheatmap(matlab_compare_df_genes,
                    cluster_cols = T,
                    color = col_pal,
                    cellheight = 3.2,
                    cellwidth = 60
)

row_order <- ph_matlab$tree_row[["order"]]
col_order <- ph_matlab$tree_col[["order"]]

ph_mdi <- pheatmap(mdi_compare_df_genes[row_order, col_order],
                   cluster_cols = F,
                   cluster_rows = F,
                   cellheight = 3.2,
                   cellwidth = 60,
                   color = col_pal
)
