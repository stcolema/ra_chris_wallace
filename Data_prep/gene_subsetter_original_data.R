#!/usr/bin/env Rscript

# Small script showing how to extract the relevant genes for the gene sets from
# Chris

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
      type = "character", default = ".csv",
      help = "file extension of target files (only used if --all set to TRUE) [default= %default]",
      metavar = "character"
    ),

    # Directory to write to
    optparse::make_option(c("-w", "--write_dir"),
      type = "character", default = ".",
      help = "directory to write files to [default= %default]",
      metavar = "character"
    ),

    # gene set to use (current options based on data from Chris)
    optparse::make_option(c("-g", "--gene_set"),
      type = "character",
      default = "all",
      help = "gene set to extract from the expression data [default= %default]", # One of `small', `medium', `big' or `all'",
      metavar = "character"
    ),

    # File to load gene data from (curently has to be .RData file)
    optparse::make_option(c("--gene_data"),
      type = "character",
      default = "stephen-genesets.RData",
      help = "data file of gene sets [default= %default]",
      metavar = "character"
    ),

    # .csv file containing probe-gene key
    optparse::make_option(c("-p", "--probe_key"),
      type = "character",
      default = "probe_key.csv",
      help = "file containing conversion key from Probe ID to Gene [default= %default]",
      metavar = "character"
    ),

    # Instruction to time programme
    optparse::make_option(c("-t", "--time"),
      type = "logical",
      default = FALSE,
      help = "instruction to print time programme takes to run [default= %default]",
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

# Not very clever - should allow user to input "all" or name
# create_gene_set <- function(cmd) {
#   if (cmd == "small") {
#     return("small")
#   }
#   if (cmd == "medium") {
#     return("med")
#   }
#   if (cmd == "big") {
#     return("big")
#   }
#   c("small", "med", "big")
# }

create_gene_set <- function(cmd, data) {
  if (cmd == "all") {
    return(names(data))
  }
  cmd
}

last_element <- function(x) {
  last_elem <- length(x)
  x[[last_elem]]
}

write_data <- function(file_name,
                       extension,
                       write_dir,
                       gene_data,
                       probe_key,
                       gene_set) {
  # Strip the files (do here as can utilise vectorised functions)
  # files_to_write <- strsplit(paste0("/", file_name), "/([^/]*).") %>%
  #   sapply(., last_element) %>%
  #   tools::file_path_sans_ext()

  files_to_write <- file_name %>%
    tools::file_path_sans_ext() %>%
    basename()

  # print(files_to_write)

  # List to record number of people and probes lost due to NA cleaning
  num_files <- length(file_name)

  # Load the dene sets data
  load(gene_data, gene_sets <- new.env())

  # Read in the probe key
  probe_key <- fread(probe_key, header = T)

  networks <- create_gene_set(gene_set, gene_sets)

  n_sets <- length(networks)
  
  # Create the sub directories to write to
  write_sub_dir <- paste0(write_dir, networks, "/")
  
  for(sub_dir in write_sub_dir){
    dir.create(sub_dir, showWarnings = FALSE)
  }
  
  # Iterate over the files - use index as can then access the read files and
  # write files both
  for (i in 1:num_files) {
    read_file <- file_name[[i]]
    write_file <- files_to_write[[i]]

    # Read in the data with the first row as the header
    dt <- fread(read_file, header = T)

    # Hold the people ID columns in their own object
    people_id_cols <- dt[, c(1, 2)]
    dt_rel <- dt[, -c(1, 2)]

    # Remove the X from the Probe IDs in the column names (for comparison with
    # the probe key file) (this is not necessary with data.table's fread)
    probes <- stringr::str_remove(colnames(dt[, -c(1, 2)]), "X")

    for(i in 1:n_sets){
      
      set <- networks[[i]]
      curr_dir <- write_sub_dir[[i]]
      
    # for (set in networks) {
      # curr_write_file <- paste0(write_file, "_", set)

      # Extract the relevant subset
      dt_subset <- dt[dt$gene_name %in% gene_sets[[set]]$external_gene_name, ]

      # Find the column indices of the gene network
      network_indices <- which(probe_key$Gene %in% gene_sets[[set]]$external_gene_name)
      
      # network_indices <- na.omit(match(
      #   gene_sets[[set]]$external_gene_name,
      #   probe_key$Gene
      # ))

      # Find which of these probes is present in the dataset
      cols_to_keep <- na.omit(match(probe_key$ProbeID[network_indices], probes))

      # Select this subset
      gene_set_dt <- dt_rel[, ..cols_to_keep]

      # Bind the output with the people id columns
      dt_out <- cbind(people_id_cols, gene_set_dt)

      # Write to a csv file
      curr_write_file <- paste0(curr_dir, "/", write_file, "_", set, ".csv")
      data.table::fwrite(dt_out, file = curr_write_file)
    }
  }
}

# Read in arguments from the command line
args <- input_arguments()

# Initial time
stm_i <- Sys.time()

# Read in expression data files
files <- read_in_data(args)

# Create the directroy if it's does not already exist
dir.create(args$write_dir, showWarnings = FALSE)

# Write files
write_data(
  files, args$extension, args$write_dir,
  args$gene_data,
  args$probe_key,
  args$gene_set
)

stm_o <- Sys.time()

if (args$time) {
  print(stm_o - stm_i)
}

# This is where I ran the specific case the above script is based upon
if (F) {
  load("Data/stephen-genesets.RData", gene_sets <- new.env())
  ls.str(gene_sets)

  # Not really necessary
  big_gene_set <- gene_sets$bignet
  mid_gene_set <- gene_sets$midnet
  small_gene_set <- gene_sets$smallnet

  # Read in example data - this should be done after transofrmation
  my_data <- read.csv("/home/MINTS/sdc56/Desktop/MDI/Data/VSN_NA_data/CD4.csv")
  my_data_2 <- fread("/home/MINTS/sdc56/Desktop/MDI/Data/Original_data/CD4_GE_Corrected4_Covars.txt", header = T, sep = " ")

  #
  probes <- stringr::str_remove(colnames(my_data_2[, -c(1, 2)]), "X")

  # Read in the probe key
  # probe_key <- read.csv("Analysis/probe_key.csv")
  probe_key <- read.csv("Data/full_probe_key.csv")

  my_data_3 <- my_data_2[, -c(1, 2)]
  ncol(my_data_3)


  # which(small_gene_set$external_gene_name %in% probe_key$Gene)
  
  network_indices <- which(probe_key$Gene %in% small_gene_set$external_gene_name)
  
  # network_indices <- na.omit(match(small_gene_set$external_gene_name, probe_key$Gene))
  cols_to_keep <- na.omit(match(probe_key$ProbeID[network_indices], probes))

  small_set_dt <- my_data_3[, ..cols_to_keep]

  out_small <- cbind(my_data_2[, 1:2], small_set_dt)

  # Add in this new variable
  my_data$gene_name <- probe_key$Gene[match(my_data$V1, probe_key$ProbeID)]

  # Extract the relevant subset
  rel_data <- my_data[my_data$gene_name %in% small_gene_set$external_gene_name, ]
}
