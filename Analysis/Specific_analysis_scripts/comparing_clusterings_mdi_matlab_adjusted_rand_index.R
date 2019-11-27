

library(tibble)
library(ggplot2)
library(mclust)
library(magrittr)
library(pheatmap)

file_path <- "/home/MINTS/sdc56/Desktop/MDI_small_geneset_outputs/matlab_output_no_vsn_no_norm/"

# Tibbles of interest
many_seeds_5000_tbl_filename <- paste0(
  "/home/MINTS/sdc56/Desktop/MDI_small_geneset_outputs/",
  "Many_seeds_small_7_output_5000/compare_tibble.rds"
)

matlab_no_norm_445_tbl_filename <- paste0(
  "/home/MINTS/sdc56/Desktop/MDI_small_geneset_outputs/",
  "matlab_output_no_vsn_no_norm/compare_tibble.rds"
)


probe_key_file <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/probe_key.csv"
probe_key <- data.table::fread(probe_key_file)

many_seeds_5000_tbl <- readRDS(many_seeds_5000_tbl_filename)
matlab_no_norm_445_tbl <- readRDS(matlab_no_norm_445_tbl_filename)

# The datasets
datasets <- many_seeds_5000_tbl$dataset

# Various numbers for sizing and looping
n_datasets <- length(datasets)
n_seeds <- nrow(many_seeds_5000_tbl$mdi_allocation[[1]])
n_comparisons <- n_datasets * n_seeds

# Create the vector of dataset names
dataset_vec <- lapply(datasets, rep, n_seeds) %>% unlist()

# Create the vector of seeds
seed_vec <- rep(1:n_seeds, n_datasets)

# Data frame to hold adjusted rand index for a given seed of C++ MDI and a given
# dataset compared to the MATLAB clustering for the same dataset
rand_df <- data.frame(
  Dataset = dataset_vec,
  Seed = seed_vec,
  Rand_index = rep(0, n_comparisons)
)
# Iterate and make rand index
for (i in 1:n_comparisons) {
  
  # The current dataset
  curr_dataset <- rand_df$Dataset[i]
  
  # Current seed for C++
  curr_seed <- rand_df$Seed[i]

  # The entry in the MATLAB tibble
  matlab_tibble_index <- matlab_no_norm_445_tbl$dataset == curr_dataset

  # The MATLAB clustering
  curr_matlab_clustering <- matlab_no_norm_445_tbl$pred_allocation %>%
    magrittr::extract(matlab_tibble_index) %>%
    magrittr::extract2(1)

  # The Many seeds tibble row entry
  mason_tibble_index_outer <- many_seeds_5000_tbl$dataset == curr_dataset
  
  # Within the MDI allocations the row to compare to the MATLAB
  mason_tibble_index_inner <- curr_seed

  # Extract the MDI custering for the approrpiate dataset and seed
  curr_mason_clustering <- many_seeds_5000_tbl$mdi_allocation %>%
    magrittr::extract(mason_tibble_index_outer) %>%
    magrittr::extract2(1) %>%
    magrittr::extract(mason_tibble_index_inner, ) %>%
    unlist()

  # Compute the addjusted rand index
  curr_rand_index <- adjustedRandIndex(curr_matlab_clustering, curr_mason_clustering)

  # Record
  rand_df$Rand_index[i] <- curr_rand_index
}

# Plot
ggplot(data = rand_df, aes(x = Rand_index)) +
  geom_histogram() +
  facet_wrap(~Dataset) +
  labs(
    title = "Distribution of adjusted rand index between C++ and MATLAB MDI",
    subtitle = "Comparing final clustering for 5,000 different seeds of C++ MDI to MATLAB",
    x = "Adjusted rand index",
    y = "Count"
  )
