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
main_wd <- "~/Desktop/Transposed_data_na_0.1/Gene_subsets/"
setwd(main_wd)
plt_type <- "pdf" # one of "pdf", "jpeg", "png",


sub_dir <- c("Big_set", "Med_set", "Small_set")

gene_sets <- c("big", "med", "small")

gene_set_dir <- paste0(main_wd, sub_dir)

num_sets <- length(sub_dir)

files_present <- vector("list", num_sets)
names(files_present) <- gene_sets

for(i in 1:num_sets){
  files_present[[i]] <- list.files(path = gene_set_dir[[i]],  pattern=".csv")
}

dataset_names <- strsplit(files_present[[1]], "_") %>% lapply(., `[[`, 1) %>% unlist()

num_files <- files_present %>% lapply(length)
total_num_files <- num_files %>% unname() %>% unlist() %>% sum()

data_tibble <- tibble(Gene_set = rep(names(num_files), unlist(num_files)),
  Tissue = rep(dataset_names, num_sets)
)

do_pca <- T

genes_present <- c()

# data_sub_lst <- vector("list", )
data_lst <- vector("list", num_sets)
names(data_lst) <- gene_sets
# data_lst %>% lapply(., `<-`, list())

for(set_name in gene_sets){
  data_lst[[set_name]] <- vector("list", num_files[[set_name]])
}


data_lst <- vector("list", total_num_files)

# %>%
#   lapply(., `<-`, list())
# 
# `<-`(x1, 3)


# Put all the data in a list of data tables
inner_ind <- 0
n_files <- 0
for(i in 1 :num_sets){
  
  inner_ind <- inner_ind + n_files
  
  curr_set <- gene_sets[[i]]
  curr_dir <- gene_set_dir[[i]]
  n_files <- length(files_present[[i]])
  for(j in 1 : n_files){
  # for(f in files_present[[i]]){
    
    # print(curr_set, j)
    f <- files_present[[i]][[j]]
    
    data_set <- f
    
    data_lst[[inner_ind + j]] <- fread(paste(curr_dir, f, sep = "/"), header = T)
    # data_lst[[curr_set]][[j]] <- fread(paste(curr_dir, f, sep = "/"), header = T)
    # data_lst[[f]] <- fread(f, header = T)
    # data_tibble[i, ]$Data <- fread(paste(curr_dir, f, sep = "/"), header = T)
  }
}

# data_lst %>% as.data.frame() %>% str()

data_tibble$Data <- data_lst

# for (f in file_name) {
#   data_lst[[f]] <- fread(f, header = T)
# }

# Acquire the relevant file names
files_to_write <- files_present %>% lapply(., tools::file_path_sans_ext)

# files_to_write <- strsplit(names(data_lst), "_?_(.*?)_?") %>% 
#   lapply("[[", 2) %>%
#   unlist()

num_datasets <- length(data_lst)


pca_lst <- list()
pca_plot_lst <- list()

pca_plt_dir <- paste(gene_set_dir, "PCA_plots", sep = "/")

inner_ind_pca <- 0
n_files_pca <- 0

# plt_type <- "png"

for (i in 1:num_sets) {
  
  curr_files_to_write <- files_to_write[[i]]
  curr_files_to_read <- files_present[[i]]
  
  curr_plt_dir <- pca_plt_dir[[i]]
  
  inner_ind_pca <- inner_ind_pca + n_files_pca
  
  n_files_pca <- length(curr_files_to_read)
  for(j in 1:n_files_pca){
  
# Carry out PCA and record the biplot
    if (do_pca) {
      raw_file_name <- file_name[[j]]
      edit_file_name <- curr_files_to_write[[j]]
      
      pca_lst[[inner_ind_pca + j]] <- .res_pca <- prcomp(t(data_tibble[inner_ind_pca + j,]$Data[[1]][, -1])) # exclude gene_name
      
      # pca_lst[[edit_file_name]] <- .res_pca <- prcomp(t(data_lst[[raw_file_name]][, -1]))
      
      pca_title <- paste0(edit_file_name, ": PCA for individuals")
      
      pca_plot_lst[[inner_ind_pca + j]] <- fviz_pca_ind(.res_pca,
        col.ind = "contrib", 
        # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
        title = pca_title
      ) +
        scale_color_gradient2(
          low = "black", mid = "blue",
          high = "red", midpoint = 1.0
        )
      
      save_name <- paste(curr_plt_dir, pca_title, sep = "/")
      ggsave(save_name, pca_plot_lst[[inner_ind_pca + j]], device = plt_type)
    }
  }
}

# data_tibble$Data[1][[1]]


summary(data_tibble$Data[13])

data_tibble$PCA <- pca_lst
data_tibble$PCA_plot <- pca_plot_lst

cd15_data <- data_tibble %>% 
  dplyr::filter(Tissue == "CD15")

# Express number of 0's present as percentage
colSums(cd15_data$Data[[1]] ==0)/nrow(cd15_data$Data[[1]] )*100


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
