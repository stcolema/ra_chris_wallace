
library(optparse)
library(magrittr)
library(data.table)

# User inputs from command line
input_arguments <- function() {
  option_list <- list(
    # Directory to read from
    optparse::make_option(c("-d", "--dir"),
                          type = "character",
                          default = ".",
                          help = "directory to read files from [default= %default]",
                          metavar = "character"
    ),
    # Directory to read from
    optparse::make_option(c("-w", "--write_dir"),
                          type = "character",
                          default = NULL,
                          help = "directory to write files to [default= %default]",
                          metavar = "character"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

# Read in arguments
args <- input_arguments()
data_dir <- args$dir
write_dir <- args$write_dir

# If NULL values for write directory, write to same location as read directory
if(is.null(write_dir)){
  write_dir <- data_dir
}

# If it does not already exist, create the write directory
dir.create(write_dir, showWarnings = FALSE)

# List the files of interest 
curr_names <- list.files(path = data_dir, full.names = T, include.dirs = F) %>%
  grep("csv", ., value = TRUE)

# Possible file names (lacking the extension)
possible_names <- c(
  "CD14",
  "CD15",
  "CD4",
  "CD8",
  "CD19",
  "IL",
  "PLA",
  "RE",
  "TR"
  ) 

# The output file names
new_names <- possible_names %>% paste0(write_dir, "/", ., ".csv")

# Numbers for FOR loops
n_files <- length(curr_names)
n_options <- length(possible_names)

# Loop over files writing a copy to the new location and name
for(i in 1:n_files){
  
  curr_name <- curr_names[[i]]
  curr_file <- fread(curr_name)
  
  # Find the new name that matches the current file
  for(j in 1 : n_options){
    
    n <- possible_names[[j]]
    
    # If the current new name is contained within the current file, find the 
    # corresponding new save path
    if(grepl(n, curr_name)){
      out_name <- new_names[[j]]
    }
  }
  
  # Write the new file
  fwrite(curr_file, out_name)
}
