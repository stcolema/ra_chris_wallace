#!/usr/bin/env Rscript

# Script to generate several new datasets based on old ones bt with more diffuse
# clusters

# === Libraries ================================================================

function_dir <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/Analysis_script_functions/"

function_scripts <- c(
  "plot_similarity_matrices.R",
  "plot_comparison_expression_clustering.R",
  "create_psm.R"
)

for (f in paste0(function_dir, function_scripts)) {
  source(f)
}

# For tibbles
library(tibble) # for dataframe of lists

# For data wrangling
library(dplyr)

# Heatmapping
library(pheatmap) # install.packages("pheatmap", dep = T)

# Colour palettes
library(RColorBrewer)

# For MVN
library(MASS)

# for fwrite
library(data.table)

# For command line inputs
library(optparse)

# === Functions ================================================================

input_arguments <- function() {
  option_list <- list(

    # tibble to use
    optparse::make_option(c("-t", "--tibble"),
      type = "character",
      default = NULL,
      help = "File name of tibble to extract cluster and expression data from.",
      metavar = "character"
    ),

    # Directory to write to
    optparse::make_option(c("-d", "--dir"),
      type = "character",
      default = ".",
      help = "Directory to write files to [default= %default]",
      metavar = "character"
    ),

    # Instruction to create directory to write to
    optparse::make_option(c("-c", "--create_dir"),
      type = "logical",
      default = TRUE,
      help = "Instruction to create directory to write to [default= %default]",
      metavar = "logical"
    ),

    # Instruction to create directory to write to
    optparse::make_option(c("-n", "--n_gen"),
      type = "integer",
      default = 10L,
      help = "Number of data point sto generate [default= %default]",
      metavar = "integer"
    ),


    # Range of weights to use for generating data
    optparse::make_option(c("-w", "--weight_range"),
      type = "double",
      default = c(0:5 / 5),
      help = "Range of weights to use for generating data [default= %default]",
      metavar = "double"
    ),

    # Label of first cluster to use
    optparse::make_option(c("--cluster_label_1"),
      type = "integer",
      default = NA,
      help = "Label of first cluster to use [default= %default]",
      metavar = "integer"
    ),

    # Label of second cluster to use
    optparse::make_option(c("--cluster_label_2"),
      type = "integer",
      default = NA,
      help = "Label of second cluster to use [default= %default]",
      metavar = "integer"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

find_sub_cluster_indices <- function(all_names, names_1, names_2) {

  # Find the relevant indices for the two clusters
  cl_1_ind <- match(names_1, all_names)
  cl_2_ind <- match(names_2, all_names)
  cl_new_ind <- c(cl_1_ind, cl_2_ind)

  list(
    "cl_1_ind" = cl_1_ind,
    "cl_2_ind" = cl_2_ind,
    "cl_new_ind" = cl_new_ind
  )
}

find_sub_clusters <- function(cl_alloc, expr_data, cl_label_1, cl_label_2) {
  all_names <- names(cl_alloc)

  names_1 <- all_names[which(cl_alloc == cl_label_1)]
  names_2 <- all_names[which(cl_alloc == cl_label_2)]

  sub_ind <- find_sub_cluster_indices(all_names, names_1, names_2)

  # Collect the expression data for each cluster
  cl_1 <- expr_data[sub_ind$cl_1_ind, ]
  cl_2 <- expr_data[sub_ind$cl_2_ind, ]
  cl_new <- expr_data[sub_ind$cl_new_ind, ]

  list(
    "cl_1" = cl_1,
    "cl_2" = cl_2,
    "cl_new" = cl_new
  )
}

generate_data <- function(cl_alloc, expr_data, cl_label_1, cl_label_2,
                          w_1 = 1,
                          w_2 = 1,
                          n_gen = 20,
                          new_names = NULL) {
  sub_cl <- find_sub_clusters(cl_alloc, expr_data, cl_label_1, cl_label_2)

  gen_data <- generate_clusters(sub_cl$cl_1, sub_cl$cl_2, sub_cl$cl_new,
    w_1 = w_1,
    w_2 = w_2,
    n_gen = n_gen,
    new_names = new_names
  )
}

generate_clusters <- function(cl_1, cl_2, cl_new,
                              w_1 = 1,
                              w_2 = 1,
                              n_gen = 20,
                              new_names = NULL) {
  if (n_gen < 1) {
    stop("`n_gen` argument must be strictly positive.")
  }

  # Find the cluster means and covariances
  mu_1 <- apply(cl_1, 2, mean)
  mu_2 <- apply(cl_2, 2, mean)

  # Covariances
  # cov_1 <- cov(cl_1)
  # cov_2 <- cov(cl_2)

  # The cluster sizes
  n_1 <- nrow(cl_1)
  n_2 <- nrow(cl_2)
  n <- n_1 + n_2

  # Define the mean to sample with - allow it to be moved towards one cluster or
  # the other
  mu_new <- (w_1 * n_1 * mu_1 + w_2 * n_2 * mu_2) / (w_1 * n_1 + w_2 * n_2)

  # Define the new clusters covariance
  cov_new <- cov(cl_new)

  # Weighted covariance for new cluster
  weight_vec <- c(rep(w_1, n_1), rep(w_2, n_2))
  cov_new_lst <- cov.wt(cl_new, wt = weight_vec)

  cov_new <- cov_new_lst$cov
  mu_new_2 <- cov_new_lst$center

  if (any(mu_new_2 - mu_new > 1e-6)) {
    stop("Calculated mean by me not equal to mean from `cov.wt` function.")
  }

  # Generate mew data (n = 1 causes issues)
  if (n_gen == 1) {
    if (is.null(new_names)) {
      new_names <- "New_gene_1"
    }

    gen_data <- mvrnorm(n_gen, mu_new, cov_new) %>%
      as.data.frame() %>%
      t() %>%
      set_rownames(new_names)
  } else {
    if (is.null(new_names)) {
      new_names <- paste0("New_gene_", 1:n_gen)
    }

    gen_data <- mvrnorm(n_gen, mu_new, cov_new) %>%
      set_rownames(new_names)
  }

  gen_data
}

create_multiple_datasets <- function(expr_data_list,
                                     cl_alloc_list,
                                     weight_range,
                                     n_datasets,
                                     cl_label_1,
                                     cl_label_2,
                                     n_gen = 10,
                                     names_1 = NULL,
                                     names_2 = NULL,
                                     dataset_names = NULL,
                                     save_dir = "./Generated_data/") {
  if (n_gen %% 2 != 0) {
    stop("n_gen must be divisble by 2.")
  }

  if (n_gen < 2) {
    stop("n_gen must be strictly positive.")
  }

  if (is.null(dataset_names)) {
    dataset_names <- paste0("D", 1:n_datasets)
  }

  dir.create(save_dir, showWarnings = FALSE)
  gen_save_dir <- paste0(save_dir, "weights_")

  # The number to generate for each cluster
  n_gen_ind <- n_gen / 2

  # If no names given generate some
  if (is.null(names_1)) {
    names_1 <- paste0("New_genes_", 1:n_gen_ind)
  }

  if (is.null(names_2)) {
    names_2 <- paste0("New_genes_", (n_gen_ind + 1):(2 * n_gen_ind))
  }

  # Loop over weights
  for (w in weight_range) {
    loc_dir <- paste0(gen_save_dir, w, "/")
    dir.create(loc_dir, showWarnings = FALSE)

    # Loop over datasets
    for (i in 1:n_datasets) {
      curr_dataset <- dataset_names[i]

      f_name <- paste0(loc_dir, curr_dataset, ".csv")

      expr_data <- expr_data_list[[i]]
      cl_alloc <- cl_alloc_list[[i]]

      # Generate data weighted towards cluster 1 (if w > 1, else towards cluster 2)
      new_data_1 <- generate_data(cl_alloc, expr_data, cl_label_1, cl_label_2,
        w_1 = w,
        w_2 = 1,
        n_gen = n_gen_ind,
        new_names = names_1
      )

      # Generate data weighted towards cluster 2 (if w > 1, else towards cluster 1)
      new_data_2 <- generate_data(cl_alloc, expr_data, cl_label_1, cl_label_2,
        w_1 = 1,
        w_2 = w,
        n_gen = n_gen_ind,
        new_names = names_2
      )

      # Combine the new data
      gen_data <- rbind(new_data_1, new_data_2)

      new_expr_dataset <- rbind(expr_data, gen_data)

      # Write the new data to .csv
      write.csv(new_expr_dataset, f_name, row.names = TRUE)
    }
  }
}

# === Main script ==============================================================

# Read in arguments
args <- input_arguments()

# The filename of the tibble to use
my_tibble_file <- args$tibble

if (is.null(my_tibble_file)) {
  stop("Please set a valid file path for `--tibble` input.")
}

# Read in the tibble
my_tibble <- readRDS(my_tibble_file)

# Create the save directory if instructed to
write_dir <- args$dir
create_dir <- args$create_dir
if (create_dir) {
  dir.create(write_dir, showWarnings = FALSE)
}

# Inputs defining clusters to generate
n_gen <- args$n_gen
weight_range <- args$weight_range

# Find the dataset names and the number thereof from the tibble
dataset_names <- unique(my_tibble$dataset)
n_datasets <- length(dataset_names)

# Labels of clusters to generate around
cl_label_1 <- args$cluster_label_1
cl_label_2 <- args$cluster_label_2

# Create the datasets
create_multiple_datasets(my_tibble$expression_data,
  my_tibble$pred_allocation,
  weight_range,
  n_datasets,
  cl_label_1,
  cl_label_2,
  n_gen = n_gen,
  names_1 = NULL,
  names_2 = NULL,
  dataset_names = dataset_names,
  save_dir = write_dir
)
