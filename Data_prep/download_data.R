#!/usr/bin/env Rscript

# Rscript to download and save data from CEDAR cohort

library(data.table)
library(R.utils)
library(magrittr)
library(optparse)

# User inputs from command line
input_arguments <- function() {
  option_list <- list(
    
    # Directroy to write files to
    optparse::make_option(c("-d", "--dir"),
                          type = "character",
                          default = "./",
                          help = "directory to write files to [default= %default]", 
                          metavar = "character"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

# Read in arguments
args <- input_arguments()

save_dir <- args$dir

# Create the directroy if it's does not already exist
dir.create(save_dir, showWarnings = FALSE)

# Website for probe expression data for CEDAR cohort
data_url <- paste0(
  "http://139.165.108.18/srv/genmol/permanent/",
  "1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/",
  "CEDAR_GE/GE_Corr/"
)

# Individual files
files <- c("TR_GE_Corrected4_Covars.txt.gz",
  "CD14_GE_Corrected4_Covars.txt.gz",
  "CD19_GE_Corrected4_Covars.txt.gz",
  "CD4_GE_Corrected4_Covars.txt.gz",
  "IL_GE_Corrected4_Covars.txt.gz",
  "CD8_GE_Corrected4_Covars.txt.gz",
  "CD15_GE_Corrected4_Covars.txt.gz",
  "RE_GE_Corrected4_Covars.txt.gz",
  "PLA_GE_Corrected4_Covars.txt.gz"
)

# Combine files with url for full read string
read_files <- paste0(data_url, files)

# Create files to save download to
write_file <- paste0(save_dir, files)

# Number of files to iterate over
num_files <- length(files)

# Download and save files
for(i in 1: num_files){
  download.file(read_files[[i]], write_file[[i]])
}

# data <- read.csv(
#   gzfile(write_file[[i]]),
#   sep= " ",
#   header=TRUE,
#   stringsAsFactors=FALSE)
# names(data)[1] <- sub("X\\.","",names(data)[1])
# 
# head(data[, 1:4])
