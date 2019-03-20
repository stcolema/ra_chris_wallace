#!/usr/bin/env Rscript

# Example of a call:
# Rscript ra_chris_wallace/Data_prep/matrix_transposer.R -a TRUE
# -d ./MDI/Data/Original\ data/ -w ./VSN_NA_min_data/ -n TRUE
# --na_people_threshold 0.1 --na_probe_threshold 0.1 -t TRUE -v TRUE


# Rscript to apply variance stabilisation

# For command line arguments
library(optparse) # install.packages("optparse")

# Load data.table to access fread and fwrite functions
library(data.table) # install.packages("data.table", dep = T)

# Load magrittr for the pipe %>%
library(magrittr)

# Load dplyr for access to select function
# library(dplyr)

# For variance stabilisation functions
# library(vsn, quietly = T, verbose = F)

# User inputs from command line
input_arguments <- function() {
  option_list <- list(

    # File to convert (if all is TRUE this is not used)
    optparse::make_option(c("-f", "--file"),
      type = "character",
      default = NA,
      help = "dataset file name",
      metavar = "character"
    ),

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("-a", "--all"),
      type = "logical",
      default = FALSE,
      help = "command to transpose all [EXT] files in currect directory [default= %default]",
      metavar = "logical"
    ),

    # Directory to read from
    optparse::make_option(c("-d", "--dir"),
      type = "character",
      default = ".",
      help = "directory to read files from (used if all set to TRUE) [default= %default]",
      metavar = "character"
    ),

    # File extension to be accepting
    optparse::make_option(c("-e", "--extension"),
      type = "character",
      default = ".csv",
      help = "file extension of target files (only used if --all set to TRUE) [default= %default]",
      metavar = "character"
    ),

    # Directory to write to
    optparse::make_option(c("-w", "--write_dir"),
      type = "character",
      default = ".",
      help = "directory to write files to [default= %default]",
      metavar = "character"
    ),

    # Instruction to time programme
    optparse::make_option(c("-t", "--time"),
      type = "logical", default = FALSE,
      help = "instruciton to record runtime of function [default= %default]",
      metavar = "logical"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

read_in_data <- function(args) {
  # If instructed to process all files, find all suitable files (i.e. ending in file_ext)
  if (args$all) {
    files_present <- list.files(path = args$dir)
    file_name <- grep(args$extension, files_present, value = TRUE) %>%
      paste0(args$dir, "/", .)
    return(file_name)
  }
  # If not transforming all files check if the user supplied a filename and
  # check if the file exists (I think this should not be working - not sure why
  # as I do not tell it where to look)

  supplied_file <- paste0(args$dir, "/", args$file)
  if (!is.na(args$file) & file.exists(supplied_file)) {
    return(supplied_file)
  } else {
    if (!is.na(args$file) & !file.exists(supplied_file)) {
      cat("File not found. \n")
      cat(paste(
        "Please check that the read directory (--dir) and filename",
        "(--file) are correct. \n"
      ))
    }
    stop("Script aborted.")
  }
  file_name
}

# Function to apply variance stabilisation to data
vsn_data_table <- function(dt, exponent = 2) {

  # Remove the id column, remove any empty observations
  # Convert to matrix, return to base values (as log2)
  # Apply vsn
  dt %>%
    .[rowSums(.) != 0, ] %>%    # Remove any empty observations
    as.matrix() %>%             # Convert to matrix (required by vsnMatrix)
    exponent^. %>%              # Exponentiate to move from log scale
    vsn::vsnMatrix() %>%        # Apply variance stabilization
    vsn::exprs() %>%            # get the expression data
    data.table::data.table()    # put it in a data table
}

# Function to apply to file names to extract relevant part
last_element <- function(x) {
  last_elem <- length(x)
  x[[last_elem]]
}


write_data <- function(file_name, extension, write_dir) {

  # Strip the files (do here as can utilise vectorised functions)
  files_to_write <- basename(tools::file_path_sans_ext(file_name))

  # List to record number of people and probes lost due to NA cleaning
  num_files <- length(file_name)
  dim_drop <- data.table::data.table(
    File = files_to_write,
    N_probes_lost = rep(0, num_files),
    N_people_lost = rep(0, num_files)
  )

  # Iterate over the files - use index as can then access the read files and
  # write files both
  for (i in 1:num_files) {
    read_file <- file_name[[i]]
    write_file <- files_to_write[[i]]

    # Read in the data with the first row as the header
    dt <- data.table::fread(read_file, header = T)


    # Variance stabilisation
    keep_names <- names(dt) != "V1"

    # don't normalise the probe ids - remove them (note data.table's odd format)
    vsn_dt <- vsn_data_table(dt[, ..keep_names])

    # Add back in the probe ids
    vsn_dt$V1 <- dt$V1

    # Rearragne order with V1 (the probe ids) in the first position
    col_order <- names(vsn_dt) %>%
      .[. != "V1"] %>%
      c("V1", .)

    dt_out <- data.table::setcolorder(vsn_dt, col_order)


    # Write to a csv file
    write_file <- paste0(write_dir, "/vsn_", write_file, ".csv")
    data.table::fwrite(dt_out, file = write_file)
  }
}

# Call functions to receive arguments and write a csv of the ransposed data
args <- input_arguments()
stm_i <- Sys.time()
files <- read_in_data(args)
write_data(files, args$extension, args$write_dir)

if (args$time) {
  print(Sys.time() - stm_i)
}
