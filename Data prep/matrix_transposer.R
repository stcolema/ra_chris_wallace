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

    # Directory to read from
    make_option(c("-d", "--dir"),
      type = "character", default = ".",
      help = "directory to read files from (used if all set to TRUE)",
      metavar = "character"
    ),

    # File extension to be accepting
    make_option(c("-e", "--extension"),
      type = "character", default = ".txt",
      help = "file extension of target files (only used if --all set to TRUE)",
      metavar = "character"
    ),

    # Directory to write to
    make_option(c("-w", "--write_dir"),
      type = "character", default = ".",
      help = "directory to write files to",
      metavar = "character"
    ),
    
    # Directory to write to
    make_option(c("-n", "--na"),
                type = "logical", default = FALSE,
                help = "instruction to remove NAs from data",
                metavar = "logical"
    ) 

    # # Index of columns containing row names
    # make_option(c("-r", "--row_names"),
    #   type = "complex", default = 0,
    #   help = "Index of columns containing row names",
    #   metavar = "complex"
    # )
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
    # cat("Enter a filename: ")
    # file_name <- readLines("stdin", 1)
    #
    # while (!file.exists(file_name)) {
    #   cat("File not found, please try again or quit [:q]: ")
    #   file_name <- readLines("stdin", 1)
    #   if (file_name == ":q") {
    #     stop("User terminated script.")
    #   }
    # }
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

write_data <- function(file_name, extension, write_dir, row_names, 
                       remove.na = FALSE) {
  for (file in file_name) {
    # Read in the data with the first row as the header
    dt <- fread(file, header = T)

    # Remove the first two columns (the ID and its duplicate) and transpose the data
    col_classes <- sapply(dt, class)
    cols_to_drop <- col_classes != "numeric"
    names_to_drop <- names(cols_to_drop) [unname(cols_to_drop)]
    conv_dt <- dt %>% dplyr::select(-c(names_to_drop)) %>%
      as.matrix() %>%
      t()
    
    if(remove.na){
      conv_dt <- conv_dt %>% na.omit()
    }
    
    # if (row_names != 0) {
    #   # Ensure that the first two columns are indeed duplicates
    #   if (any(dt[, 1] != dt[, 2])) {
    #     cat("ID columns not duplicates. Please inspect data.\n")
    #   }
    # 
    #   conv_dt <- dt[, -(row_names)]
    # }
    # conv_dt <- dt %>%
    #   as.matrix() %>%
    #   t()

    # Create list of column names
    id <- dt[, 1] %>%
      create_col_names_vector()

    # Assign the ID as column names
    colnames(conv_dt) <- id

    # Convert to data.table to use fwrite
    dt_out <- as.data.table(conv_dt)
    row.names(dt_out) <- row.names(conv_dt)

    file_name <- strsplit(paste0("/", file) , "/([^/]*).")[[1]]
    file_name <- file_name[[length(file_name)]]

    # Write to a csv file
    file_to_write <- sub(paste0(extension, "$"), "", file_name) %>%
      paste0(write_dir, "/transposed_", ., ".csv")

    fwrite(dt_out, file = file_to_write, row.names = T)
  }
}

# Call functions to receive arguments and write a csv of the ransposed data
args <- input_arguments()
files <- read_in_data(args)
write_data(files, args$extension, args$write_dir, args$row_names, args$na)