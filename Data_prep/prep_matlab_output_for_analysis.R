#!/usr/bin/env Rscript

# This file adds column names to the output of Paul Kirk's MATLAB MDI code
# corresponding to the column names for the output of Sam Mason's MDI code.
# 
# Example call:
# Rscript prep_matlab_output_for_analysis.R -o example_mdi_analysis.csv -d 
# example_input_data.csv -n 3 -f data_ready_for_analysis.csv

library(mdiHelpR)
library(data.table)
library(stringr)

input_arguments <- function() {
  option_list <- list(

    # File name of MATLAB MDI output to add colnames to
    optparse::make_option(c("--mdi_output", "-o"),
      type = "character",
      default = NA,
      help = "MATLAB MDI output file to add colnames to [default= %default]",
      metavar = "character"
    ),

    # Input data to MDI that contains the sample names
    optparse::make_option(c("--data", "-d"),
      type = "character",
      default = NA,
      help = "Input data with rownames corresponding to sample names [default= %default]",
      metavar = "character"
    ),

    # Name to save file to
    optparse::make_option(c("--file_name", "-f"),
      type = "character",
      default = "./mdi_output.csv",
      help = "File name to save to [default= %default]",
      metavar = "character"
    ),

    # The number of datasets used in this analysis
    optparse::make_option(c("--n_datasets", "-n"),
      type = "integer",
      default = NA,
      help = "Number of datasets in analysis [default= %default]",
      metavar = "character"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

args <- input_arguments()

orig_data_file <- args$data
mdi_output_file <- args$mdi_output
file_name <- args$file_name
n_datasets <- args$n_datasets

# timecourse <- read.csv("~/Desktop/MDI_data/timecourse.csv", row.names = 1)

# example_out <- read.csv("~/Desktop/Gran_unnorm/Consensus_clustering_more_clusters/gran_unnorm_consensus_clustering_large_c.csv")
# colnames(example_out)

# genes <- row.names(timecourse)
orig_data <- read.csv(orig_data_file, row.names = 1)
sample_names <- row.names(orig_data)

# n_datasets <- 3
dataset_names <- paste0("Dataset", 1:n_datasets)

# mdi_output <- read.csv("~/Desktop/ppi_harbison_marina_Granovskaia_time_course_551_101.csv",
# header = F)

mdi_output <- read.csv(mdi_output_file, header = F)

# mdi_output
# colnames(mdi_output)
# dim(mdi_output)



new_colnames <- c()
for (i in 1:n_datasets) {
  new_colnames <- c(new_colnames, paste(dataset_names[i], sample_names, sep = "_"))
}


mdi_param_names <- mdiHelpR::createParameterNames(n_datasets)
new_colnames <- c(mdi_param_names, new_colnames) %>% str_replace("-", ".")

colnames(mdi_output) <- new_colnames

# write.csv(mdi_output, "~/Desktop/matlab_output_100.csv")
write.csv(mdi_output, file_name, row.names = F)

# file_name
# any(new_colnames != colnames(example_out))
#
# mismatch <- which(new_colnames != colnames(example_out))
#
# new_colnames[1:10]
# colnames(example_out)[1:10]
#
#
# new_colnames[mismatch]
# colnames(example_out)[mismatch]
#
# ?str_replace()
#
