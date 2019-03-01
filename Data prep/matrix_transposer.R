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
    optparse::make_option(c("-f", "--file"),
      type = "character", default = NA,
      help = "dataset file name", metavar = "character"
    ),

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("-a", "--all"),
      type = "logical", default = FALSE,
      help = "command to transpose all [EXT] files in currect directory [default= %default]",
      metavar = "logical"
    ),

    # Directory to read from
    optparse::make_option(c("-d", "--dir"),
      type = "character", default = ".",
      help = "directory to read files from (used if all set to TRUE) [default= %default]",
      metavar = "character"
    ),

    # File extension to be accepting
    optparse::make_option(c("-e", "--extension"),
      type = "character", default = ".txt",
      help = "file extension of target files (only used if --all set to TRUE) [default= %default]",
      metavar = "character"
    ),

    # Directory to write to
    optparse::make_option(c("-w", "--write_dir"),
      type = "character", default = ".",
      help = "directory to write files to [default= %default]",
      metavar = "character"
    ),

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

  # Call an error if the dir input is not 1 or 2
  if(! dir %in% c(1,2)) stop("dir must be one of 1 (for rows) or 2 (for columns")
  
  # For any remaining NAs replace with the column / row median as instructed
  if (dir == 2) {
    # Remove columns / rows with too high a proportion of NAs
    dt_na_dropped <- na_threshold(dt,
      threshold = col_threshold,
      cols = TRUE
    )

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

  # If dir == 1 do the analog of the above to the rows.
  dt_na_dropped <- na_threshold(dt,
    threshold = row_threshold,
    cols = FALSE
  )

  # If we want to do this by row, it seems easiest (as this is unlikely to be of
  # actual interest) t ouse the preceding method. This means we must transpose
  # our data
  dt_no_na <- purrr::pmap(dt_na_dropped, function(x) {
    y <- ifelse(is.na(x), median(x, na.rm = TRUE), x)
  })

  dummy_dt <- dt_na_dropped %>%
    as.matrix() %>%
    t() %>%
    as.data.table()

  # As above, cycle through each collumn replacing NAs with the median value
  for (col in colnames(dummy_dt)) {
    e <- substitute(
      X := ifelse(is.na(X), median(X, na.rm = TRUE), X),
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

write_data <- function(file_name, extension, write_dir, row_names,
                       remove.na = FALSE,
                       na_people_threshold = 0.0,
                       na_probe_threshold = 0.0,
                       dir = 2) {
  for (file in file_name) {
    # Read in the data with the first row as the header
    dt <- fread(file, header = T)

    # Remove the first two columns (the ID and its duplicate) and transpose the data
    col_classes <- sapply(dt, class)
    cols_to_drop <- col_classes != "numeric"
    names_to_drop <- names(cols_to_drop) [unname(cols_to_drop)]

    # Drop the ID variables and transpose the data
    conv_dt <- dt %>%
      dplyr::select(-c(names_to_drop)) %>%
      as.matrix() %>%
      t() %>%
      as.data.table()

    # Assign row.names
    # row.names(conv_dt) <- colnames(dt[, -(1:2)])

    # Create list of column names
    id <- dt[, 1] %>%
      create_col_names_vector()

    # Assign the ID as column names
    colnames(conv_dt) <- id

    # Get the probe IDs into our data table
    conv_dt$V1 <- colnames(dt[, -(1:2)])

    # If removing NAs, do so
    if (remove.na) {
      conv_dt <- conv_dt %>%
        na_handling(
          row_threshold = na_probe_threshold,
          col_threshold = na_people_threshold,
          dir = dir
        )
    }

    # print(summary(conv_dt))

    # Convert to data.table to use fwrite
    dt_out <- as.data.table(conv_dt)
    # row.names(dt_out) <- conv_dt$V1
    # row.names(dt_out) <- row.names(conv_dt)

    # Rearragne order with V1 (the probe ids) in the first position
    col_order <- names(dt_out) %>%
      .[. != "V1"] %>%
      c("V1", .)

    dt_out <- setcolorder(dt_out, col_order)

    # Create file name
    file_name <- strsplit(paste0("/", file), "/([^/]*).")[[1]]
    file_name <- file_name[[length(file_name)]]

    # Write to a csv file
    file_to_write <- sub(paste0(extension, "$"), "", file_name) %>%
      paste0(write_dir, "/transposed_", ., ".csv")

    data.table::fwrite(dt_out, file = file_to_write)
  }
}

# Call functions to receive arguments and write a csv of the ransposed data
args <- input_arguments()
stm_i <- Sys.time()
files <- read_in_data(args)
write_data(files, args$extension, args$write_dir, args$row_names,
  remove.na = args$na,
  na_people_threshold = args$na_people_threshold,
  na_probe_threshold = args$na_probe_threshold
)

print(Sys.time() - stm_i)
