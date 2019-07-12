#!/usr/bin/env Rscript

# Example of a call:
# Rscript ra_chris_wallace/Data_prep/matrix_transposer.R -a TRUE
# -d ./MDI/Data/Original\ data/ -w ./VSN_NA_min_data/ -n TRUE
# --na_people_threshold 0.1 --na_probe_threshold 0.1 -t TRUE -v TRUE


# Rscript to convert data from CEDAR cohorts (found here: http://cedar-web.giga.ulg.ac.be/)
# Call from the command line with input arguments
# Writes to current directory

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
      default = ".txt",
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

    # # Instruction to transpose data
    # optparse::make_option(c("-t", "--transpose"),
    #   type = "logical",
    #   default = TRUE,
    #   help = "instruciton to transpose data [default= %default]",
    #   metavar = "logical"
    # ),

    # Instruction to remove NAs
    optparse::make_option(c("-n", "--na"),
      type = "logical", default = FALSE,
      help = "instruction to remove NAs from data [default= %default]",
      metavar = "logical"
    ),

    # Threshold at which to remove PROBES
    optparse::make_option(c("--na_probe_threshold"),
      type = "double", default = 0.0,
      help = "na threshold at which to remove PROBES from data [default= %default]",
      metavar = "double"
    ),

    # Threshold at which to remove PEOPLE
    optparse::make_option(c("--na_people_threshold"),
      type = "double", default = 0.0,
      help = "na threshold at which to remove PEOPLE from data [default= %default]",
      metavar = "double"
    ),

    # Instruction to apply variance stabilisation to the data
    optparse::make_option(c("-v", "--vsn"),
      type = "logical", default = FALSE,
      help = "instruction to apply variance stabilisation to the data [default= %default]",
      metavar = "logical"
    ),

    # Instruction to apply variance stabilisation to the data
    optparse::make_option(c("-c", "--centre"),
      type = "logical", default = FALSE,
      help = "instruction to centre the data [default= %default]",
      metavar = "logical"
    ),

    # Instruction to apply variance stabilisation to the data
    optparse::make_option(c("-s", "--scale"),
      type = "logical", default = FALSE,
      help = "instruction to scale the data [default= %default]",
      metavar = "logical"
    ),

    # Instruction to time programme
    optparse::make_option(c("-t", "--time"),
      type = "logical", default = FALSE,
      help = "instruciton to record runtime of function [default= %default]",
      metavar = "logical"
    )

    # # Index of columns containing row names
    # make_option(c("-r", "--row_names"),
    #   type = "complex", default = 0,
    #   help = "Index of columns containing row names",
    #   metavar = "complex"
    # )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

read_in_data <- function(args) {
  # If instructed to process all files, find all suitable files (i.e. ending in file_ext)
  if (args$all) {
    files_present <- list.files(path = args$dir) #, pattern = paste0("*", args$extension))
    
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


# Remove columns / rows with too high a proportion of NAs
na_threshold <- function(dt, threshold = 0.1, cols = TRUE) {
  if (cols) {
    cols_to_keep <- !(colSums(is.na(dt)) / nrow(dt)) > threshold

    # If this happens (that we keep all columns) data.table seems to get
    # confused and move to a 0x0 object. Thus we return dt if not removing
    # any columns
    if (sum(cols_to_keep) == ncol(dt)) {
      return(dt)
    }
    new_dt <- dt[, cols_to_keep, with = F]

    return(new_dt)
  }

  # Find all rows not exceeding our threshold for proportion of NAs
  rows_to_keep <- !(rowSums(is.na(dt)) / ncol(dt)) > threshold

  # Again ,if keeping all rows return the original object
  if (sum(rows_to_keep) == nrow(dt)) {
    return(dt)
  }
  new_dt <- dt[rows_to_keep, ]

  new_dt
}

# Convert NAs in a vector to the median value
na_to_median <- function(x) {

  # first convert each column into numeric if it is from factor
  # x <- as.numeric(as.character(x))

  # convert the item with NA to median value from the column
  x[is.na(x)] <- median(x, na.rm = TRUE)
}

# Samll function to call is using apply or purrr::map to replace NAs with median
na_replace <- function(x) {
  y <- ifelse(is.na(x), median(x, na.rm = TRUE), x)
  y
}

# Handle NAs by removing rows / columns with too high a proportion of NAs and
# replace remaining NAs with the associated column's median value
na_handling <- function(dt, row_threshold = 0.1, col_threshold = 0.1, dir = 2) {
  # dt is the data in data.table format
  # row_threshold is the threshold for which rows with proportions of NAs
  # exceeding this are removed
  # col_threshold is the corresponding vlaue for the columns
  # dir is the dimension to fill NAs with median vlaues (1 is rows, 2 columns)

  # Call an error if dir outside acceptable values
  if (!dir %in% c(1, 2)) stop("dir must be one of 1 (rows) or 2 (columns).")

  # Remove columns / rows with too high a proportion of NAs
  dt_cols_dropped <- na_threshold(dt,
    threshold = col_threshold,
    cols = TRUE
  )

  dt_na_dropped <- na_threshold(dt_cols_dropped,
    threshold = row_threshold,
    cols = FALSE
  )

  # For any remaining NAs replace with the column / row median as instructed
  if (dir == 2) {
    for (col in colnames(dt_na_dropped)) {

      # Create the expression to evaluate within the data.table
      # Specifically, replace any NAs with the column median otherwise leave
      # them untouched
      e <- substitute(
        X := ifelse(is.na(X), median(X, na.rm = TRUE), X),
        list(X = as.symbol(col))
      )

      # Evaluate this within the data.table
      dt_no_na <- dt_na_dropped[, eval(e)]
    }

    # consdier purrr (tried this - significantly slower, roughly 4 times as long)
    # dt_no_na <- purrr::map_df(dt_na_dropped, na_replace)

    return(dt_no_na)
  }

  # If we want to do this by row, it seems easiest (as this is unlikely to be of
  # actual interest) t ouse the preceding method. This means we must transpose
  # our data
  # dt_no_na <- purrr::pmap(dt_na_dropped, function(x) {
  #   y <- ifelse(is.na(x), median(x, na.rm = TRUE), x)
  # })

  dummy_dt <- dt_na_dropped %>%
    as.matrix() %>%
    t() %>%
    as.data.table()

  # As above, cycle through each collumn replacing NAs with the median value
  for (col in colnames(dummy_dt)) {
    e <- substitute(
      X := ifelse(is.na(X), min(X, na.rm = TRUE), X),
      list(X = as.symbol(col))
    )

    # Evaluate this within the data.table
    dt_no_na <- dummy_dt[, eval(e)]
  }

  # Return to the correct orientation and object
  dt_no_na %<>%
    as.matrix() %>%
    t() %>%
    as.data.table()

  # Re-assign column names
  colnames(dt_no_na) <- colnames(dt_na_dropped)

  dt_no_na
}


# Function to apply to file names to extract relevant part
last_element <- function(x) {
  last_elem <- length(x)
  x[[last_elem]]
}

# Function to apply variance stabilisation to data
vsn_data_table <- function(dt, exponent = 2) {

  # Remove the id column, remove any empty observations
  # Convert to matrix, return to base values (as log2)
  # Apply vsn
  dt %>%
    .[rowSums(.) != 0, ] %>% # Remove any empty observations
    as.matrix() %>% # Convert to matrix (required by vsnMatrix)
    # lumi::inverseVST() %>% # remove variance stabilisation transform (this is not possible with the information availabe)
    # exponent^. %>% # Exponentiate to move from log scale
    vsn::vsnMatrix() %>% # Apply variance stabilization
    vsn::exprs() %>% # get the expression data
    data.table::data.table() # put it in a data table
}


write_data <- function(file_name, extension, write_dir,
                       remove.na = FALSE,
                       na_people_threshold = 0.0,
                       na_probe_threshold = 0.0,
                       dir = 2,
                       do_vsn = FALSE,
                       scale = FALSE,
                       centre = FALSE) {

  # Strip the files (do here as can utilise vectorised functions)
  files_to_write <- strsplit(paste0("/", file_name), "/([^/]*).") %>%
    sapply(., last_element) %>%
    tools::file_path_sans_ext()

  # List to record number of people and probes lost due to NA cleaning
  num_files <- length(file_name)
  dim_drop <- data.table::data.table(
    File = files_to_write,
    N_probes_lost = rep(0, num_files),
    N_people_lost = rep(0, num_files)
  )

  print(file_name)
  
  # Iterate over the files - use index as can then access the read files and
  # write files both
  for (i in 1:num_files) {
    read_file <- file_name[[i]]
    write_file <- files_to_write[[i]]

    # Read in the data with the first row as the header
    dt <- data.table::fread(read_file, header = T)

    # Remove the first two columns (the ID and its duplicate) and transpose the data
    col_classes <- sapply(dt, class)
    cols_to_drop <- col_classes != "numeric"

    names_to_drop <- names(cols_to_drop) [unname(cols_to_drop)]

    # Drop the ID variables and transpose the data
    conv_dt <- dt %>%
      dplyr::select(-c(names_to_drop)) %>%
      as.matrix() %>%
      t() %>%
      data.table::as.data.table()

    # Assign row.names
    # row.names(conv_dt) <- colnames(dt[, -(1:2)])

    # Create list of column names
    id <- dt[, 1] %>%
      create_col_names_vector()

    # Assign the ID as column names
    colnames(conv_dt) <- id

    # Get the probe IDs into our data table
    conv_dt$V1 <- colnames(dt[, -(1:2)])

    # Initial dimensionality before cleaning
    n_col_i <- ncol(conv_dt)
    n_row_i <- nrow(conv_dt)

    # If removing NAs, do so
    if (remove.na) {
      conv_dt <- conv_dt %>%
        na_handling(
          row_threshold = na_probe_threshold,
          col_threshold = na_people_threshold,
          dir = dir
        )
    }
    
    # Dimensionality post-cleaning
    n_col_o <- ncol(conv_dt)
    n_row_o <- nrow(conv_dt)

    # Record the number of people and probes dropped for this file
    dim_drop[dim_drop$File == write_file]$N_probes_lost <- n_row_i - n_row_o
    dim_drop[dim_drop$File == write_file]$N_people_lost <- n_col_i - n_col_o

    # Convert to data.table to use fwrite
    dt_out <- data.table::as.data.table(conv_dt)

    if(scale | centre){

      # Record the current probes actually present
      final_probes_present <- dt_out$V1

      # Scale and centre if instructed
      dt_scaled <- dt_out %>% 
        extract(, V1:=NULL) %>% # Remove non-numerical column (data.table grammar)
        t() %>% # Transpose as wish to standardise genes
        scale(scale = scale, center = centre) %>% 
        t() %>% # Return to orientation of interest
        as.data.table() # Revert to a data.table
      
      # Add back in the probes
      dt_scaled$V1 <- final_probes_present
      
      # Return to the same object name as if the IF statement is passed
      dt_out <- dt_scaled
    }
    
    if (do_vsn) {
      
      cat("\nDoing variance stabilisation normalisation.\n")
      
      # Variance stabilisation
      keep_names <- names(dt_out) != "V1"

      # don't normalise the probe ids - remove them (note data.table's odd format)
      vsn_dt <- vsn_data_table(dt_out[, ..keep_names])

      # Add back in the probe ids
      vsn_dt$V1 <- dt_out$V1

      # move back to an object not unique to this if statement
      dt_out <- vsn_dt
    }

    # Rearragne order with V1 (the probe ids) in the first position
    col_order <- names(dt_out) %>%
      .[. != "V1"] %>%
      c("V1", .)

    # Set the column order and scale that data
    dt_out <- data.table::setcolorder(dt_out, col_order)

    # Write to a csv file
    write_file <- paste0(write_dir, "/transposed_", write_file, ".csv")
    data.table::fwrite(dt_out, file = write_file)
  }

  if (remove.na) {
    print(dim_drop)
    data.table::fwrite(dim_drop,
      file = paste0(write_dir, "Observations_lost.csv")
    )
  }
}

# Call functions to receive arguments and write a csv of the ransposed data
args <- input_arguments()
stm_i <- Sys.time()
files <- read_in_data(args)

# Create the directroy if it's does not already exist
dir.create(args$write_dir, showWarnings = FALSE)

write_data(files, args$extension, args$write_dir,
  remove.na = args$na,
  na_people_threshold = args$na_people_threshold,
  na_probe_threshold = args$na_probe_threshold,
  do_vsn = args$vsn,
  scale = args$scale,
  centre = args$centre
)

if (args$time) {
  print(Sys.time() - stm_i)
}
