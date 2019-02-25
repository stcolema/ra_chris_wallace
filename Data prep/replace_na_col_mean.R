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

# Load dplyr for access to select function
# library(dplyr) 

# User inputs from command line
input_arguments <- function() {
  option_list <- list(
    
    # File to convert (if all is TRUE this is not used)
    make_option(c("-f", "--file"),
                type = "character", default = NA,
                help = "dataset file name", metavar = "character"
    ),
    
    # Convert all files in target destination (default is FALSE)
    make_option(c("-a", "--all"),
                type = "logical", default = FALSE,
                help = "command to transpose all .txt files in currect directory [default= %default]",
                metavar = "logical"
    ),
    
    # File extension to be accepting
    make_option(c("-e", "--extension"),
                type = "character", default = ".csv",
                help = "file extension of target files (only used if --all set to TRUE)",
                metavar = "character"
    ),
    
    # Directory to read from
    make_option(c("-d", "--dir"),
                type = "character", default = ".",
                help = "directory to read files from (used if all set to TRUE)",
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
    file_name <- grep(args$extension, files_present, value = TRUE) %>%
      paste0(args$dir, "/", .)
    return(file_name)
  }
  # If not transforming all files check if the user supplied a filename and
  # check if the file exists (I think this should not be working - not sure why
  # as I do not tell it where to look)
  
  supplied_file <- paste0(args$dir, "/", args$file)
  if (!is.na(args$file) & file.exists(supplied_file)) {
    # file_name <- args$file
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

write_data <- function(file_name) {
  for (file in file_name) {
    # Read in the data with the first row as the header
    dt <- fread(file, header = T)
    
    col_mean <- colMeans(dt)
    na_columns <- names(dt)[colSums(is.na(dt)) != 0]
    
    dt_out <- dt[, (na_columns) := lapply(na_columns, function(x) {
      x <- get(x)
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      x
    }), by = V1]
    
    # for(i in 2:ncol(dt)){
    # for(na_col in na_columns){
    #   dt[is.na(dt[ ,..i]), ..i] <- unname(col_mean[i])
    # }

    file_name <- strsplit(paste0("/", file) , "/([^/]*).")[[1]]
    file_name <- file_name[[length(file_name)]]
    
    # Write to a csv file
    file_to_write <- file_name %>%
      paste0("na_filled_", .)
    
    fwrite(dt_out, file = file_to_write, row.names = F)
  }
}

# Call functions to receive arguments and write a csv of the ransposed data
args <- input_arguments()
files <- read_in_data(args)
write_data(files)