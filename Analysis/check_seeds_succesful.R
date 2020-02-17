#!/usr/env/bin/Rscript

# Check how many chains failed to run and failed to run the correct number of 
# iterations

# Load libraries
library(magrittr)   # chain and variosu associated functions
library(optparse)   # command line input
library(parallel)   # parallel computing

fileLength <- function(f, header = FALSE) {
  x <- read.csv(f, header = header)
  nrow(x)
}

input_arguments <- function() {
  option_list <- list(

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("-d", "--dir"),
      type = "character",
      default = "./",
      help = "directory to scan for CI output [default= %default]",
      metavar = "character"
    ),

    optparse::make_option(c("-n", "--n_seeds"),
      type = "integer",
      default = 1000,
      help = "the number of seeds one expected to have run successfully [default= %default]",
      metavar = "integer"
    ),

    optparse::make_option(c("--n_iter"),
      type = "integer",
      default = 2,
      help = "the number of iterations one expected in each .csv file [default= %default]",
      metavar = "integer"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

args <- input_arguments()
main_data_dir <- args$dir
expected_runs <- args$n_seeds
expected_samples <- args$n_iter

# main_data_dir <- "~/Documents/PhD/Year_1/Consensus_clustering/Analysis/MDI_test_data/Difficult_convergence/p_250/MATLAB_CI/"
# expected_runs <- 5000

# List subdirectories
sub_dirs <- list.dirs(main_data_dir, recursive = F)
sub_dirs_short <- list.dirs(main_data_dir, recursive = F, full.names = F)

# Number of subdirectories
n_dirs <- length(sub_dirs)

# The number of output files acutally expected
expected_files <- paste0("out_seed_", 1:expected_runs, ".csv")

# The object to hold the names of the files to rerun
seeds_to_run <- vector("list", n_dirs) %>% set_names(sub_dirs_short)

# Iterate over each subdirectory seeing a) which files are entirely missing and 
# b) which files have too few samples
for (i in 1:n_dirs) {
  curr_files <- list.files(sub_dirs[i], pattern = ".csv", full.names = F)
  missing_seeds <- expected_files[which(!(expected_files %in% curr_files))]

  my_data_length <- mclapply(paste0(sub_dirs[i], "/", curr_files), fileLength) %>%
    set_names(curr_files)

  seeds_to_rerun <- names(my_data_length)[which(unlist(my_data_length) < expected_samples)]

  seeds_to_run[[i]] <- c(missing_seeds, seeds_to_rerun)
}

# See which files are missing
print(seeds_to_run)

