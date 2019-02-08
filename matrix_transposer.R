#!/usr/bin/env Rscript

# Rscript to convert data from CEDAR cohorts (found here: http://cedar-web.giga.ulg.ac.be/)
# Call from the command line with input arguments
# Writes to current directory

# For command line arguments
library("optparse") # install.packages("optparse")

# Load data.table to access fread and fwrite functions
library(data.table) # install.packages("data.table", dep = T)

# Load magrittr for the pipe %>%
library(magrittr)

# User inputs from command line
input_arguments <- function() {
  option_list <- list(
    make_option(c("-f", "--file"),
      type = "character", default = NA,
      help = "dataset file name", metavar = "character"
    ),
    make_option(c("-a", "--all"),
      type = "logical", default = FALSE,
      help = "command to transpose all .txt files in currect directory [default= %default]",
      metavar = "logical"
    ),
    make_option(c("-d", "--dir"),
      type = "character", default = ".",
      help = "directory to read files from (used if all set to TRUE)",
      metavar = "character"
    ),
    make_option(c("-e", "--extension"),
      type = "character", default = ".txt",
      help = "file extension of target files (only used if --all set to TRUE)",
      metavar = "character"
    ),
    make_option(c("-w", "--write_dir"),
                type = "character", default = ".",
                help = "directory to write files to",
                metavar = "character"
    )
  )
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
}

read_in_data <- function(args) {
  # If instructed to process all files, find all suitable files (i.e. ending in file_ext)
  if (args$all) {
    files_present <- list.files(path = args$dir)
    file_name <- grep(args$extension, files_present, value = TRUE)
    return(file_name)
  }
  if (!is.na(args$file) & file.exists(args$file)) {
    file_name <- args$file
    return(file_name)
  } else {
    if (!is.na(args$file) & !file.exists(args$file)) {
      print("File not found.")
    }
    cat("Enter a filename: ")
    file_name <- readLines("stdin", 1)

    # file_name <- readline(prompt="Enter a filename: ")
    while (!file.exists(file_name)) {
      cat("File not found, please try again or quit [:q]: ")
      file_name <- readLines("stdin", 1)
      # quit <- readline(prompt="File not found, please try again or quit [:q]: ")
      if (file_name == ":q") {
        stop("User terminated script.")
      }
      # file_name <- readline(prompt="Enter a filename: ")
    }
  }
  file_name
}

# Function to convert a data frame to (not used as type conversion hidden)
transpose_data <- function(dt) {
  dt %>%
    as.matrix() %>%
    t()
}

# Given a n x 1 column vector in data.table type, returns a list of the entries
create_col_names_vector <- function(dt) {
  # Create column names from ID column
  dt %>%
    unname() %>%
    unlist()
}

write_data <- function(file_name, extension, write_dir) {
  for (file in file_name) {
    # Read in the data with the first row as the header
    dt <- fread("~/Desktop/Data/CD14_GE_Corrected4_Covars.txt", header = T)

    # Remove the first two columns (the ID and its duplicate) and transpose the data
    conv_dt <- dt[, -(1:2)] %>%
      as.matrix() %>%
      t()

    # ID column
    # Ensure that the first two columns are indeed duplicates
    if (any(dt[, 1] != dt[, 2])) {
      print("ID columns not duplicates. Please inspect data.")
    }

    # Create list of column names
    id <- dt[, 1] %>%
      create_col_names_vector()

    # Assign the ID as column names
    colnames(conv_dt) <- id

    # Convert to data.table to use fwrite
    dt_out <- as.data.table(conv_dt)
    row.names(dt_out) <- row.names(conv_dt)

    # Write to a csv file
    file_to_write <- sub(paste0(extension, "$"), "", file) %>% 
      paste0(write_dir, "/transposed_", ., ".csv")
    
    fwrite(dt_out, file = file_to_write, row.names = T)
  }
}

# Call functions to receive arguments and write a csv of the ransposed data
args <- input_arguments()
files <- read_in_data(args)
write_data(files, args$extension, args$write_dir)
