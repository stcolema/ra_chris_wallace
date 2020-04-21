#!/usr/env/bin/Rscript

################################################################################
#                                                                              #
# Title: Consensus inference performance                                       #
#                                                                              #
# Aim: Analyse saved clustering output of Mason's implementation of MDI for    #
# Consensus inference simulations run April 2020.                              #
# Compare models to the ground truth using ARI and a normalised Frobenius      #
# product.                                                                     #
#                                                                              #
# Author: Stephen Coleman                                                      #
# Date: 21/04/2020                                                             #
#                                                                              #
################################################################################

# For my ggplot theme settings
library(mdiHelpR)

# File name interactions
library(stringr)

# Plotting
library(gglot2)

# I don't think this is used
library(tibble)

# To read in data
library(data.table)

# Visualising heatmaps (no longer used)
library(pheatmap)

# For frobenius.prod
library(matrixcalc)

# For arandi and maxpear for model evaluation and predicted clustering
library(mcclust)

# For facet_wrap_paginate
library(ggforce)

# For the saving of Sparse matrices
library(Matrix)

# For command line arguments
library(optparse)

# === Functions ================================================================

# The Frobenius product is used to measure model uncertainty
normalisedFrobeniusProduct <- function(sm, truth) {
  frobenius.prod(sm, truth) / frobenius.prod(truth, truth)
}

inputArguments <- function() {
  option_list <- list(

    # Data to cluster
    optparse::make_option(c("-d", "--dir"),
      type = "character",
      default = "./",
      help = "File path to model output [default= %default]",
      metavar = "character"
    ),

    # Save directory
    optparse::make_option(c("-s", "--save_dir"),
      type = "character",
      default = "./",
      help = "Directory to save model and predicted clustering to [default= %default]",
      metavar = "character"
    ),

    # Expected chain length
    optparse::make_option(c("-l", "--length"),
      type = "integer",
      default = 10001,
      help = "Expected chain length; will call an error if any chains are less than this [default= %default]",
      metavar = "integer"
    ),

    # Thinning factor applied within each chain
    optparse::make_option(c("-t", "--thin"),
      type = "character",
      default = "1 10 100 1000",
      help = "Thinning factor applied within each chain [default= %default]",
      metavar = "character"
    ),

    # Burn-in to apply to each chain
    optparse::make_option(c("--n_seeds"),
      type = "integer",
      default = 100,
      help = "Number of chains in analysis [default= %default]",
      metavar = "integer"
    ),

    # Burn-in to apply to each chain
    optparse::make_option(c("--sim_num"),
      type = "integer",
      default = 1,
      help = "Current simulation [default= %default]",
      metavar = "integer"
    ),

    # Burn-in to apply to each chain
    optparse::make_option(c("--truth_dir"),
      type = "character",
      default = NA,
      help = "Directory containing truth [default= %default]",
      metavar = "character"
    ),

    # Columns of interest for convergence
    optparse::make_option(c("-c", "--cols"),
      type = "character",
      default = "-1",
      help = "Columns to subset by in ``fread`` [default= %default]",
      metavar = "character"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

# === Initialisation ===========================================================

args <- inputArguments()

# Directory to read from and to save Sparse matrices to
data_dir <- args$dir
save_dir <- args$save_dir
truth_dir <- args$truth_dir

dir.create(save_dir)

# Details of analysis
expected_length <- args$length
thin <- args$thin %>% 
  str_split(" ") %>% 
  unlist() %>% 
  as.numeric()

n_seeds <- args$n_seeds
sim_num <- args$sim_num

# Columns to exclude as irrelevant to clustering
subset_cols <- args$cols %>% 
  str_split(" ") %>% 
  unlist() %>% 
  as.numeric()

# For reproducibility
set.seed(1)

# Use a subset of people in the tests
# Remember to set subset cols to -1 to exclude Mass_Parameter
subset_cols <- -1
if (interactive()) {
  subset_cols <- c(2:600)
  n_people <- length(subset_cols)
}

# Files containing simulation data and structure
cluster_id_file <- paste0(truth_dir, "cluster_IDs_", sim_num, ".Rds")
orig_data_file <- paste0(truth_dir, "dataset_", sim_num, ".csv")

# True labelling and data being clustered
truth <- readRDS(cluster_id_file)
orig_data <- read.csv(orig_data_file)

# True coclustering matrix
true_cc <- createSimilarityMat(matrix(truth[1:n_people], nrow = 1))

# List input files for reading in and for anlaying
files_full <- list.files(data_dir, full.names = T) %>%
  str_sort(numeric = T)

files_short <- list.files(data_dir, full.names = F) %>%
  str_sort(numeric = T)

# Find the details of each file based on a naming convention of:
# "${Model type}N{Number of iterations}T{thinning factor}S{random seed}.csv"
file_details <- files_short %>%
  str_extract_all("[1-9][0-9]*") %>%
  unlist() %>%
  as.numeric() %>%
  matrix(ncol = 3, byrow = T) %>%
  set_colnames(c("n", "thin", "seed")) %>%
  as.data.frame()

# Inspect the file details
# file_details %>%
#   head()
#
# file_details %>% nrow()

# Seed order (NOT USED)
seed_order <- file_details$seed %>%
  order()

# Files ordered by number of iterations
iter_order <- file_details$n %>%
  order()

# Inspect
# file_details[seed_order, ] %>%
# head()

# Number of seeds and iterations present
n_seeds <- max(file_details$seed)
n_iter <- max(file_details$n)

# Iterations saved across all consensus inference models in simulations run April 2020
iter_used <- c(
  seq(1, 10, 1),
  seq(20, 100, 10),
  seq(200, 1000, 100),
  seq(2000, 10000, 1000)
)

# Used for book-keeping
n_model <- length(iter_used)


# Inspect
# files_full[seed_order] %>%
# head()

# Find dimensions of datasets by reading in a 0-row data.table
s1 <- fread(files_full[1], select = subset_cols, nrows = 0)
n_people <- ncol(s1)

# Old object to hold psms
# n_seed_df <- n_iter_df <- matrix(0, nrow = n_people, ncol = n_people)

# seed_psm_array <- array(0, c(n_people, n_people, n_seeds, n_model))

# Number of iterations and seeds to record for model evaluation
iters_to_record <- c(1, 10, 100, 1000, 10000)
seeds_to_record <- c(1, 10, 30, 50, 100)

# The number of models recorded under each measure
n_iters_recorded <- length(iters_to_record)
n_seeds_recorded <- length(seeds_to_record)

# The data.frame that holds the score of each model
results_df <- data.frame(
  N_iter = rep(iters_to_record, n_seeds_recorded),
  N_seeds = rep(seeds_to_record, each = n_iters_recorded),
  ARI = 0,
  Frobenius_product = 0
)

# NOT USED
# recorded_values <- c("N_iter", "N_seeds", "ARI", "Frobenius_product")

# The array that hold the consensus matrix for each model
seed_psm_array <- array(0, c(n_people, n_people, n_model))

# === Analysis =================================================================

for (i in 1:n_seeds) {
  for (j in 1:n_model) {

    # Current number of iterations used
    curr_iter <- iter_used[j]

    # Which consensus output matches current requirements of seed and chain
    # length
    file_used <- which(file_details$n > curr_iter & file_details$seed == i)[1]

    # What thinning is present in the current file (influences which sample is
    # kept, if in a file where we want the thousandth sample and there is
    # thinning of 100, we want the (1000/100 + 1)=11th sample, with the +1
    # correcting for Mason always saving the 0th iteration)
    thin_present <- file_details$thin[file_used]

    # Read in the relevant samples
    curr_samples <- fread(files_full[[file_used]], select = subset_cols, nrows = 1 + curr_iter / thin_present) %>%
      as.matrix(ncol = n_people)

    # Create a coclustering matrix based on the current sample
    curr_psm <- createSimilarityMat(matrix(curr_samples[1 + curr_iter / thin_present, ], nrow = 1))

    # Update the consensus matrix
    seed_psm_array[, , j] <- .curr_cm <- seed_psm_array[, , j] * (i - 1) / i + curr_psm * (1 / i)

    # If current model is to to be saved, perform inference
    if (curr_iter %in% iters_to_record & i %in% seeds_to_record) {

      # Save a sparse representation of the similarity matrix
      sparse_cm <- Matrix(.curr_cm, sparse = TRUE)
      cm_file <- paste0(save_dir, "consensusMatrixN", curr_iter, "S", i, ".txt")
      writeMM(sparse_cm, file = cm_file)

      # Predicted clustering
      curr_pred_cl <- maxpear(.curr_cm)$cl

      # Record frobenius product and ARI
      record_ind <- which(results_df$N_iter == curr_iter & results_df$N_seeds == i)
      results_df$ARI[record_ind] <- arandi(curr_pred_cl, truth[1:n_people])
      results_df$Frobenius_product[record_ind] <- normalisedFrobeniusProduct(.curr_cm, true_cc)
    }
  }
}

# === Plotting =================================================================

write.csv(results_df, file = paste0(save_dir, "ResultsCI.csv"))

# Set labels for facet wrapping
iter_labels <- c(paste0("Number of iterations: ", results_df$N_iter))
names(iter_labels) <- results_df$N_iter

seed_labels <- c(paste0("Number of chains: ", results_df$N_seeds))
names(seed_labels) <- results_df$N_seeds

# Plots
p_ari_iter <- results_df %>%
  ggplot(aes(x = N_seeds, y = ARI)) +
  geom_line() +
  facet_wrap(~N_iter, labeller = labeller(N_iter = iter_labels), ncol = 1) # +
# labs(
#   title = "Performance stabilitses by 100 iterations",
#   subtitle = "CI for large n, large k, small distance between means",
#   x = "Number of chains used"
# )

p_ari_seed <- results_df %>%
  ggplot(aes(x = N_iter, y = ARI)) +
  geom_line() +
  facet_wrap(~N_seeds, labeller = labeller(N_seeds = seed_labels), ncol = 1) # +
# labs(
# title = "Performance stabilitses by 100 iterations",
# subtitle = "CI for large n, large k, small distance between means",
# x = "Number of iterations used"
# )


p_unc_iter <- results_df %>%
  ggplot(aes(x = N_seeds, y = Frobenius_product)) +
  geom_line() +
  facet_wrap(~N_iter, labeller = labeller(N_iter = iter_labels), ncol = 1) # +
# labs(
#   title = "Performance stabilitses after 40 iterations",
#   subtitle = "CI for large n, large k, small distance between means",
#   x = "Number of chains used"
# )

p_unc_seed <- results_df %>%
  ggplot(aes(x = N_iter, y = Frobenius_product)) +
  geom_line() +
  facet_wrap(~N_seeds, labeller = labeller(N_seeds = seed_labels), ncol = 1)

p_ari_iter + p_ari_seed
p_unc_iter + p_unc_seed

p_ari_iter + p_unc_iter

#
# results_df %>%
#   ggplot(aes(x = N_iter, y = Uncertainty)) +
#   geom_point() +
#   facet_wrap_paginate(~N_seeds, ncol = 3, nrow = 3, page =3)
#
# results_df %>%
#   ggplot(aes(x = N_seeds, y = Uncertainty)) +
#   geom_point() +
#   facet_wrap_paginate(~N_iter, ncol = 4, nrow = 4, page =1) +
# labs(
#   title = "Performance stabilitses after 40 iterations",
#   subtitle = "CI for large n, large k, small distance between means",
#   x = "Number of chains used"
# )
#
# results_df %>%
#   ggplot(aes(x = N_iter, y = Uncertainty)) +
#   geom_point() +
#   facet_wrap_paginate(~N_seeds, ncol = 4, nrow = 4, page =1) +
#   labs(
#     title = "Performance stabilitses after 40 iterations",
#     subtitle = "CI for large n, large k, small distance between means",
#     x = "Number of chains used"
#   )
#
# results_df %>%
#   ggplot(aes(x = N_seeds, y = ARI)) +
#   geom_point() +
#   facet_wrap_paginate(~N_iter, ncol = 4, nrow = 4, page =1,
#                       labeller = labeller(N_iter = facet_labels)) +
#   labs(
#     title = "Performance stabilitses after 40 iterations",
#     subtitle = "CI performance for large n, large k, small distance between means (ARI)",
#     x = "Number of chains used"
#   )
#
#
# results_df %>%
#   ggplot(aes(x = N_seeds, y = Uncertainty)) +
#   geom_point() +
#   facet_wrap_paginate(~N_iter, ncol = 4, nrow = 4, page =3,
#                       labeller = labeller(N_iter = facet_labels)) +
#   labs(
#     title = "Performance stabilitses after 40 iterations",
#     subtitle = "CI performance for large n, large k, small distance between means (Uncertainty)",
#     x = "Number of chains used"
#   )
#
#
#
# results_df %>%
#   ggplot(aes(x = N_iter, y = ARI)) +
#   geom_point() +
#   facet_wrap_paginate(~N_seeds, ncol = 4, nrow = 4, page =1) +
#   labs(
#     title = "Performance stabilitses after 40 iterations",
#     subtitle = "CI for large n, large k, small distance between means",
#     x = "Number of chains used"
#   )
#
#
#
# results_df %>%
#   ggplot(aes(x = N_seeds, y = Uncertainty)) +
#   geom_point() +
#   facet_wrap_paginate(~N_iter, ncol = 3, nrow = 3, page =2)
#
# results_df %>%
#   ggplot(aes(x = N_seeds, y = Uncertainty)) +
#   geom_point() +
#   facet_wrap_paginate(~N_iter, ncol = 3, nrow = 3, page = 3)
#
# 
# 
# 
# sim_col <- simColPal()
# my_breaks <- defineBreaks(sim_col, lb = 0)
# seed_psm_array[, , 37] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
# 
# seed_psm_array[, , 1] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
# 
# seed_psm_array[, , 10] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
# 
# seed_psm_array[, , 19] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
# 
# 
# seed_psm_array[, , 28] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
# 
# 
# my_c <- maxpear(seed_psm_array[, , 6])$cl
# arandi(my_c, truth[1:n_people])

# annotatedHeatmap(orig_data[1:n_people, 2:5], my_c)
#
# my_c <- maxpear(seed_psm_array[, , 20])$cl
# annotatedHeatmap(orig_data[1:n_people, 2:5], my_c)
# 
# seed_psm_array[, , 1] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
# 
# seed_psm_array[, , 1, 1] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
# 
# seed_psm_array[, , 10, 1] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)
# 
# 
# seed_psm_array[, , 10, 1] %>%
#   pheatmap(color = sim_col, breaks = my_breaks)