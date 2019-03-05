#!/usr/bin/env Rscript

# Load data.table to access fread and fwrite functions
library(data.table) # install.packages("data.table", dep = T)

# Load magrittr for the pipe %>%
library(magrittr)

# For select, filter
library(dplyr) # install.packages("tidyverse", dep = T)

prep_data <- function(dt) {
  row.names(dt) <- dt[, 1] %>%
    unlist()

  dt %<>%
    select(-V1)
}

# Read in data
# files_present <- list.files(path = args$dir)
# file_name <- grep(args$extension, files_present, value = TRUE) %>%
#   paste0(args$dir, "/", .)
setwd("~/Desktop/My_end/")
files_present <- list.files(path = "~/Desktop/My_end/")
file_name <- grep(".csv", files_present, value = TRUE)

data_lst <- list()
setwd("~/Desktop/My_end")

# Put all the data in a list of data tables
for(f in file_name){
  
  data_lst[[f]] <- fread(f, header = T)
  
  # Record the genes present in all datasets
  if(length(data_lst) == 1){
    genes_present  <- data.frame(V1 = data_lst[[f]]$V1)
  }
  if(length(data_lst) > 1){
    genes_present <- semi_join(
      genes_present,
      select(data_lst[[f]], V1),
      by = "V1"
    ) 
  }
}

# Move the genes present to a list
genes_present %<>% unname() %>% unlist()

# For each dataset filter by genes present in all and write to file
for(f in names(data_lst)){
  data_lst[[f]] %<>%
    filter(V1 %in% genes_present)

  file_write <- strsplit(f, "_?_(.*?)_?")[[1]][2] %>% paste0(., ".csv")
  
  fwrite(data_lst[[f]], file = file_write, row.names = F)
}