#!/usr/bin/env Rscript

# Load data.table to access fread and fwrite functions
library(data.table) # install.packages("data.table", dep = T)

# Load magrittr for the pipe %>%
library(magrittr)

# For select, filter
library(dplyr) # install.packages("tidyverse", dep = T)

# library("devtools") # install.packages("devtools", dep = T)
# install_github("kassambara/factoextra")
library(factoextra) # install.packages("factoextra", dep = T)

prep_data <- function(dt) {
  row.names(dt) <- dt[, 1] %>%
    unlist()

  dt %<>%
    select(-V1)
}

# Used in initial idea for removing points from PCA that failed to meet some 
# criterion
contrib_cut <- function(pca_res, cut = 1.5, dims = 1:3, criterion = "contrib"){
  contrib <- get_pca_ind(pca_res)[[criterion]]
  ind_to_remove <-  apply(contrib[,dims], 1, function(r) any(r > cut))
}

get_cut_data <- function(pca_lst, threshold = c(1.0, 1.5, 2.0), criterion = "contrib"){
  cut_data <- list()
  for(cut in threshold){
    cut_data[[as.character(cut)]] <- pca_lst %>% 
      lapply(., contrib_cut, cut = cut, criterion = criterion) %>% 
      lapply(., sum) %>% 
      unlist()
    # unlist(lapply(lapply(pca_lst, contrib_cut), sum))
  }
  cut_data
}


# Read in data
setwd("~/Desktop/MDI/Data/Fill_NAs_Min")
files_present <- list.files(path = "~/Desktop/MDI/Data/Fill_NAs_Min")
file_name <- grep("Covars.csv", files_present, value = TRUE)

do_pca <- T

genes_present <- c()

pca_lst <- list()
pca_plot_lst <- list()


# Put all the data in a list of data tables
for (f in file_name) {
  data_lst[[f]] <- fread(f, header = T)
}

# Acquire the relevant file names
files_to_write <- strsplit(names(data_lst), "_?_(.*?)_?") %>% 
  lapply("[[", 2) %>%
  unlist()

num_datasets <- length(data_lst)

for (i in 1:num_datasets) {
# Carry out PCA and record the biplot
  if (do_pca) {
    raw_file_name <- file_name[[i]]
    edit_file_name <- files_to_write[[i]]
    
    pca_lst[[edit_file_name]] <- .res_pca <- prcomp(t(data_lst[[raw_file_name]][, -1]))
    
    pca_title <- paste0(edit_file_name, ": PCA for individuals")
    
    pca_plot_lst[[edit_file_name]] <- fviz_pca_ind(.res_pca,
      col.ind = "contrib", 
      # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
      title = pca_title
    ) +
      scale_color_gradient2(
        low = "black", mid = "blue",
        high = "red", midpoint = 1.0
      )
    
  }
}

# If PCA was done, print the plots (notice that the first compoen)
# These are all the same - something has gone wrong
if (do_pca) {
  print(pca_plot_lst$CD14)
  print(pca_plot_lst$CD15)
  print(pca_plot_lst$CD19)
  print(pca_plot_lst$CD4)
  print(pca_plot_lst$CD8)
  print(pca_plot_lst$IL)
  print(pca_plot_lst$PLA)
  print(pca_plot_lst$RE)
  print(pca_plot_lst$TR)
}

# Dropping outliers
cleaned_data <- vector("list", 9) 
# cleaned_data <- list(length = 9)
names(cleaned_data) <- files_to_write

print(pca_plot_lst$CD14)
cd14_cols_to_drop <- c("IPC154", "IPC155")
keep_cols_cd14 <- !(colnames(data_lst$transposed_CD14_GE_Corrected4_Covars.csv) %in% cd14_cols_to_drop)
cleaned_data$CD14 <- data_lst$transposed_CD14_GE_Corrected4_Covars.csv[, ..keep_cols_cd14 ]

print(pca_plot_lst$CD15)
cd15_cols_to_drop <- c("IPC137") # c("IPC300", "IPC097", "IPC244", "IPC315", "IPC334", "IPC332", "IPC262", "IPC137")
keep_cols_cd15 <- !(colnames(data_lst$transposed_CD15_GE_Corrected4_Covars.csv) %in% cd15_cols_to_drop)
cleaned_data$CD15 <- data_lst$transposed_CD15_GE_Corrected4_Covars.csv[, ..keep_cols_cd15 ]

print(pca_plot_lst$CD19)
cd19_cols_to_drop <- c()
keep_cols_cd19 <- !(colnames(data_lst$transposed_CD19_GE_Corrected4_Covars.csv) %in% cd19_cols_to_drop)
cleaned_data$CD19 <- data_lst$transposed_CD19_GE_Corrected4_Covars.csv[, ..keep_cols_cd19 ]

print(pca_plot_lst$CD4)
cd4_cols_to_drop <- c("IPC329", "IPC331")
keep_cols_cd4 <- !(colnames(data_lst$transposed_CD4_GE_Corrected4_Covars.csv) %in% cd4_cols_to_drop)
cleaned_data$CD4 <- data_lst$transposed_CD4_GE_Corrected4_Covars.csv[, ..keep_cols_cd4 ]

print(pca_plot_lst$CD8)
cd8_cols_to_drop <- c("IPC078", "IPC049", "IPC048", "IPC050")
keep_cols_cd8 <- !(colnames(data_lst$transposed_CD8_GE_Corrected4_Covars.csv) %in% cd8_cols_to_drop)
cleaned_data$CD8 <- data_lst$transposed_CD8_GE_Corrected4_Covars.csv[, ..keep_cols_cd8 ]

print(pca_plot_lst$IL)
il_cols_to_drop <- c("IPC434")
keep_cols_il <- !(colnames(data_lst$transposed_IL_GE_Corrected4_Covars.csv) %in% il_cols_to_drop)
cleaned_data$IL <- data_lst$transposed_IL_GE_Corrected4_Covars.csv[, ..keep_cols_il ]

print(pca_plot_lst$PLA)
pla_cols_to_drop <- c("IPC029") # maybe c("IPC029", "IPC106")
keep_cols_pla <- !(colnames(data_lst$transposed_PLA_GE_Corrected4_Covars.csv) %in% pla_cols_to_drop)
cleaned_data$PLA <- data_lst$transposed_PLA_GE_Corrected4_Covars.csv[, ..keep_cols_pla ]

print(pca_plot_lst$RE)
re_cols_to_drop <- c("IPC361") # maybe c("IPC361", "IPC253")
keep_cols_re <- !(colnames(data_lst$transposed_RE_GE_Corrected4_Covars.csv) %in% re_cols_to_drop)
cleaned_data$RE <- data_lst$transposed_RE_GE_Corrected4_Covars.csv[, ..keep_cols_re ]

print(pca_plot_lst$TR)
tr_cols_to_drop <- c()
keep_cols_tr <- !(colnames(data_lst$transposed_TR_GE_Cortrcted4_Covars.csv) %in% tr_cols_to_drop)
cleaned_data$TR <- data_lst$transposed_TR_GE_Corrected4_Covars.csv[, ..keep_cols_re ]

write_dir <- "/home/MINTS/sdc56/Desktop/new_data/"

for(i in 1:num_datasets){
  file_name <- paste0(write_dir, names(cleaned_data)[[i]], ".csv")
  fwrite(cleaned_data[[i]], file = file_name, row.names = F)
}
