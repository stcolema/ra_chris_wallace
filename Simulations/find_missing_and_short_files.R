#!/usr/env/bin/Rscript

################################################################################
#                                                                              #
# Title: Find missing samples                                                  #
#                                                                              #
# Find files with missing samples or entirely absent files for Mason MDI runs. #
# Author: Stephen Coleman                                                      #
# Date: 21/04/2020                                                             #
################################################################################

# Command line
library(optparse)

# For plotting
library(ggplot2)

# For pipe, filter and mutate
suppressMessages(library(dplyr))

# For string order (unnecessary)
library(stringr)

# For fread
library(data.table)

# Pipe
library(magrittr)

# -- Functions -----------------------------------------------------------------

# Function to take arguments from the command line using the ``optparse`` library
inputArguments <- function() {
  option_list <- list(

    # Data to cluster
    optparse::make_option(c("-d", "--dir"),
      type = "character",
      default = "./",
      help = "File path to read consensus analysis from  [default= %default]",
      metavar = "character"
    ),

    # Save directory
    optparse::make_option(c("-s", "--save_dir"),
      type = "character",
      default = "./",
      help = "Directory to save plot to [default= %default]",
      metavar = "character"
    ),

    # Inference type
    optparse::make_option(c("--inference_type"),
      type = "character",
      default = "Consensus",
      help = "Inference type in question (Bayesian or Consensus) [default= %default]",
      metavar = "character"
    ),

    # Save directory
    optparse::make_option(c("-i", "--iterations_used"),
      type = "character",
      default = "10 100 1000 10001",
      help = "Number of iterations used as input in current analysis [default= %default]",
      metavar = "character"
    ),

    # Save directory
    optparse::make_option(c("-t", "--thin"),
      type = "character",
      default = "1 10 100 1000",
      help = "Thinning factors applied (pasired with iterations_used) [default= %default]",
      metavar = "character"
    ),

    # Save directory
    optparse::make_option(c("-n", "--n_seeds"),
      type = "integer",
      default = 100L,
      help = "Number of chains run [default= %default]",
      metavar = "integer"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

# === Main script ==============================================================

setMyTheme()
set.seed(1)

# Read-in command line arguments
args <- inputArguments()

# Directory to check for missing/short files
my_dir <- args$dir

# Inference performed (changes expected file names, naming convention is:
# "${Inference type}N${number of iterations}T${thinning factor}Seed${seed}.csv",
# e.g. "ConsensusN10001T1000Seed21.csv")
inference_type <- args$inference_type

# Iterations used across samples (can be string with space separating iterations 
# run); paired with thin
iterations_used <- args$iterations_used %>% 
  str_split(" ") %>%
  unlist() %>%
  as.numeric()

# Thinning factor similar to above; expected to be the same length 
thin <- args$thin %>%
  str_split(" ") %>%
  unlist() %>%
  as.numeric()

# Number of different chains run for each (n_iter, thin) pairing
n_seeds <- args$n_seeds

# Expected chain lengths
lengths_expect <- ceiling(iterations_used / thin)

# Number of thinning factors applied
n_thin <- length(thin)

# Number of iterations run applied
n_iterations_used <- length(iterations_used)

# List the files present in both long and short form
short_files <- list.files(my_dir) %>%
  stringr::str_sort(numeric = T)

long_files <- list.files(my_dir, full.names = T) %>%
  stringr::str_sort(numeric = T)

n_files <- length(long_files)

# The details of the generated files
file_details <- short_files %>%
  str_extract_all("[1-9][0-9]*") %>%
  unlist() %>%
  as.numeric() %>%
  matrix(ncol = 3, byrow = T) %>%
  set_colnames(c("n", "thin", "seed")) %>%
  as.data.frame() %>%
  mutate(
    "Expected_length" = ceiling(n / thin),
    Thin = factor(thin),
    N = factor(n)
  )

# The data.frame of the details of the expected files
expected_files <- data.frame(
  Expected_length = rep(lengths_expect, each = n_seeds),
  n = factor(rep(iterations_used, each = n_seeds)),
  thin = factor(rep(thin, each = n_seeds)),
  seed = rep(1:n_seeds, n_iterations_used)
)

# Create the vector of expected files
expected_file_names <- paste0(
  inference_type,
  "N",
  expected_files$n,
  "T",
  expected_files$thin,
  "Seed",
  expected_files$seed,
  ".csv"
)

# Check if any are missing
missing_files <- expected_file_names[which(is.na(match(expected_file_names, short_files)))]
n_files_missing <- length(missing_files)

# Find which files are shorter than expected
file_details$samples_record <- 0
for (i in 1:n_files) {
  f <- long_files[i]
  file_details$samples_record[i] <- nrow(fread(f, select = 1L))
}

# If any are missing print their name
if (n_files_missing) {
  cat(paste0("\n", n_files_missing, " files are missing. Please generate:\n"))
  cat(paste(missing_files, collapse = "\n"))
}

# If any are short print their name
n_files_short <- sum(file_details$Expected_length != file_details$samples_record)
if (n_files_short) {
  cat(paste0("\n\n", n_files_short, " files are shorter than expected. Re-run analysis for:\n"))
  f_ind <- which(file_details$Expected_length != file_details$samples_record)
  cat(paste(short_files[f_ind], collapse = "\n"))
}
