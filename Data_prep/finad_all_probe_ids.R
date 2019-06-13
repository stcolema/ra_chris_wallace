#!/usr/bin/env Rscript

library(optparse)

# Load data.table to access fread and fwrite functions
library(data.table) # install.packages("data.table", dep = T)

# Load magrittr for the pipe %>%
library(magrittr)

# For select, filter
library(dplyr) # install.packages("tidyverse", dep = T)

# library("devtools") # install.packages("devtools", dep = T)
# install_github("kassambara/factoextra")
# library(factoextra) # install.packages("factoextra", dep = T)

prep_data <- function(dt) {
  row.names(dt) <- dt[, 1] %>%
    unlist()

  dt %<>%
    select(-V1)
}

# Used in initial idea for removing points from PCA that failed to meet some
# criterion
contrib_cut <- function(pca_res, cut = 1.5, dims = 1:3, criterion = "contrib") {
  contrib <- get_pca_ind(pca_res)[[criterion]]
  ind_to_remove <- apply(contrib[, dims], 1, function(r) any(r > cut))
}

get_cut_data <- function(pca_lst, threshold = c(1.0, 1.5, 2.0), criterion = "contrib") {
  cut_data <- list()
  for (cut in threshold) {
    cut_data[[as.character(cut)]] <- pca_lst %>%
      lapply(., contrib_cut, cut = cut, criterion = criterion) %>%
      lapply(., sum) %>%
      unlist()
    # unlist(lapply(lapply(pca_lst, contrib_cut), sum))
  }
  cut_data
}

input_arguments <- function() {
  option_list <- list(

    # Directory to read from
    optparse::make_option(c("-d", "--dir"),
      type = "character",
      default = ".",
      help = "directory to read files from (used if all set to TRUE) [default= %default]",
      metavar = "character"
    ),

    # File extension to be accepting
    optparse::make_option(c("-v", "--vsn_applied"),
      type = "logical",
      default = FALSE,
      help = "logical indicating if vsn has been applied to data [default= %default]",
      metavar = "logical"
    ),

    # File extension to be accepting
    optparse::make_option(c("-m", "--use_min"),
      type = "logical",
      default = FALSE,
      help = "logical indicating if fill value should be minimum expressed value in expression data [default= %default]",
      metavar = "logical"
    ),

    # File extension to be accepting
    optparse::make_option(c("-z", "--use_zero"),
      type = "logical",
      default = FALSE,
      help = "logical indicating if fill value should be 0 [default= %default]",
      metavar = "logical"
    ),

    # File extension to be accepting
    optparse::make_option(c("-r", "--use_random"),
      type = "logical",
      default = FALSE,
      help = "logical indicating if fill value should be draws from a standard normal [default= %default]",
      metavar = "logical"
    ),

    # File extension to be accepting
    optparse::make_option(c("-o", "--order_data"),
      type = "logical",
      default = FALSE,
      help = "logical indicating if data should be ordered by Probe ID [default= %default]",
      metavar = "logical"
    ),


    # File extension to be accepting
    optparse::make_option(c("-i", "--index"),
      type = "integer",
      default = 2,
      help = "integer indicating section of string to use [default= %default]",
      metavar = "integer"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

args <- input_arguments()
dir_of_interest <- args$dir # "~/Desktop/subset_data/Small/" #
vsn_data <- args$vsn_applied # F

# Which fill to use
use_min <- args$use_min
use_zero <- args$use_zero
use_random <- args$use_random

# Order data by Probe ID
order_data <- args$order_data

if (sum(use_min, use_zero, use_random) > 1) {
  stop("Only one of 'use_min', 'use_zero' or 'use_random' may be set to TRUE.")
}

if (sum(use_min, use_zero, use_random) != 1) {
  stop("Please set one of 'use_min', 'use_zero' or 'use_random' to TRUE.")
}

name_ind <- args$index # 2 + vsn_data

# Read in data
setwd(dir_of_interest)
files_present <- list.files(path = dir_of_interest, pattern = ".csv")

# We should never be doing this anymore
if (vsn_data) {
  file_name <- grep("vsn_*", files_present, value = TRUE)
} else {
  file_name <- grep("*.csv", files_present, value = TRUE)
}

eda <- F

data_lst <- list()

genes_present <- c()

mean_lst <- list()
sd_lst <- list()

# Put all the data in a list of data tables
for (f in file_name) {
  data_lst[[f]] <- fread(f, header = T)

  if (order_data) {
    data_lst[[f]] <- data_lst[[f]] %>%
      .[order(.$V1), ]
  }


  # Record the genes present in all datasets
  genes_present <- unique(c(genes_present, data_lst[[f]]$V1))

  if (eda) {
    mean_lst[[f]] <- apply(data_lst[[f]][, -1], 2, mean)
    sd_lst[[f]] <- apply(data_lst[[f]][, -1], 2, sd)
  }
}

# Acquire the relevant file names
files_to_write <- names(data_lst) %>%
  tools::file_path_sans_ext() %>%
  basename() %>%
  paste0(., "_filled")


num_datasets <- length(data_lst)

# Move the genes present to a list
genes_present %<>% unname() %>% unlist()

empty_probes_dt <- data.table(matrix(NA,
  nrow = length(genes_present),
  ncol = length(file_name) + 1
))


colnames(empty_probes_dt) <- c("V1", files_to_write)
empty_probes_dt$V1 <- genes_present

# === Fill missing genes =======================================================

# For each dataset filter by genes present in all and write to file
for (i in 1:num_datasets) {
  f <- names(data_lst)[[i]]
  genes_to_add <- genes_present[!genes_present %in% data_lst[[f]]$V1]

  col_names <- colnames(data_lst[[f]])

  n_rows <- length(genes_to_add)
  n_cols <- length(col_names)

  # Fill the added rows with requested values
  if (use_min) {
    fill_vlaue <- min(data_lst[[f]]) - sd(as.matrix(data_lst[[f]]))
  }
  if (use_zero) {
    fill_value <- 0
  }

  if (use_min | use_zero) {
    additional_rows <- data.frame(matrix(fill_value,
      nrow = n_rows,
      ncol = n_cols
    ))
  }
  else {
    additional_rows <- matrix(rnorm(n_rows * n_cols, mean = 0, sd = 1), n_rows, n_cols) %>%
      as.data.frame()
  }

  colnames(additional_rows) <- col_names
  additional_rows$V1 <- genes_to_add

  # add the additional, empty probes and arrange in a common order
  dt_out <- data_lst[[f]] %>%
    bind_rows(additional_rows) %>%
    .[match(genes_present, .$V1), ]

  if (any(dt_out$V1 != genes_present)) {
    stop("Check order is correct")
  }



  file_write <- files_to_write[[i]] %>%
    paste0(., ".csv")
  empty_probes <- !dt_out$V1 %in% genes_to_add
  empty_probes_dt[[files_to_write[[i]]]] <- empty_probes

  # file_write %<>% paste0(., ".csv")

  fwrite(dt_out, file = file_write, row.names = F)
}

fwrite(empty_probes_dt, file = "probes_present_per_dataset.csv")

# summary(dt_out[, 1:5])

genes_present_dt <- data.table(Probe_ID = genes_present)
fwrite(genes_present_dt, "probe_IDs_present.csv")
