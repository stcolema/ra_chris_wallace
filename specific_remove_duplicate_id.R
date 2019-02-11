#!/usr/bin/env Rscript

# Load data.table to access fread and fwrite functions
library(data.table) # install.packages("data.table", dep = T)

# Load magrittr for the pipe %>%
library(magrittr)

# For select
library(dplyr) # install.packages("tidyverse", dep = T)

dt <- fread("~/Desktop/Data/CD14_GE_Corrected4_Covars.txt", header = T) %>% 
  dplyr::select(-FID)

fwrite(dt, file = "~/Desktop/My_end/CD14_GE_Corrected4_Covars_one_id.csv", row.names = T)
  