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
# files_present <- list.files(path = args$dir)
# file_name <- grep(args$extension, files_present, value = TRUE) %>%
#   paste0(args$dir, "/", .)
setwd("~/Desktop/MDI/Data/Fill_NAs_Min")
# setwd("~/Desktop/Na_filled_data/")
files_present <- list.files(path = "~/Desktop/MDI/Data/Fill_NAs_Min")
# files_present <- list.files(path = "~/Desktop/Na_filled_data/")
file_name <- grep("Covars.csv", files_present, value = TRUE)

do_pca <- T
eda <- F

# file_name <-  file_name[-7]
data_lst <- list()
# setwd("~/Desktop/My_end")

genes_present <- c()

pca_lst <- list()
pca_plot_lst <- list()

mean_lst <- list()
sd_lst <- list()


# Put all the data in a list of data tables
for (f in file_name) {
  data_lst[[f]] <- fread(f, header = T)

  # Record the genes present in all datasets
  genes_present <- unique(c(genes_present, data_lst[[f]]$V1))

  if (eda) {
    mean_lst[[f]] <- apply(data_lst[[f]][, -1], 2, mean)
    sd_lst[[f]] <- apply(data_lst[[f]][, -1], 2, sd)
  }

  # Carry out PCA and record the biplot
  # if (do_pca) {
  #   pca_lst[[f]] <- prcomp(t(data_lst[[f]][, -1]))
  # 
  #   pca_title <- paste()
  #   
  #   pca_plot_lst[[f]] <- fviz_pca_ind(pca_lst[[f]],
  #     geom.ind = "point", # show points only (nbut not "text")
  #     col.ind = "contrib"
  #   ) +
  #     scale_color_gradient2(
  #       low = "black", mid = "blue",
  #       high = "red", midpoint = 4
  #     )
  # }
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
    
    # pca_plot_lst[[edit_file_name]] <- fviz_pca_ind(pca_lst[[edit_file_name]],
    #                                   geom.ind = "point", # show points only (nbut not "text")
    #                                   col.ind = "contrib",
    #                                   title = pca_title
    # ) +
    #   scale_color_gradient2(
    #     low = "black", mid = "blue",
    #     high = "red", midpoint = 1.5
    #   )
  }
}

if (eda) {
  hist(mean_lst$transposed_CD14_GE_Corrected4_Covars.csv)
  hist(mean_lst$transposed_CD15_GE_Corrected4_Covars.csv)
  hist(mean_lst$transposed_CD19_GE_Corrected4_Covars.csv)
  hist(mean_lst$transposed_CD4_GE_Corrected4_Covars.csv)
  hist(mean_lst$transposed_CD8_GE_Corrected4_Covars.csv)
  hist(mean_lst$transposed_PLA_GE_Corrected4_Covars.csv)
  hist(mean_lst$transposed_IL_GE_Corrected4_Covars.csv)
  hist(mean_lst$transposed_RE_GE_Corrected4_Covars.csv)
  hist(mean_lst$transposed_TR_GE_Corrected4_Covars.csv)

  hist(sd_lst$transposed_CD14_GE_Corrected4_Covars.csv)
  hist(sd_lst$transposed_CD15_GE_Corrected4_Covars.csv)
  hist(sd_lst$transposed_CD19_GE_Corrected4_Covars.csv)
  hist(sd_lst$transposed_CD4_GE_Corrected4_Covars.csv)
  hist(sd_lst$transposed_CD8_GE_Corrected4_Covars.csv)
  hist(sd_lst$transposed_PLA_GE_Corrected4_Covars.csv)
  hist(sd_lst$transposed_IL_GE_Corrected4_Covars.csv)
  hist(sd_lst$transposed_RE_GE_Corrected4_Covars.csv)
  hist(sd_lst$transposed_TR_GE_Corrected4_Covars.csv)
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

# Move the genes present to a list
genes_present %<>% unname() %>% unlist()

empty_probes_dt <- data.table(matrix(NA,
                  nrow = length(genes_present),
                  ncol = length(file_name) + 1
))


colnames(empty_probes_dt) <- c("V1", files_to_write)
empty_probes_dt$V1 <-  genes_present


# For each dataset filter by genes present in all and write to file
for (f in names(data_lst)) {
  genes_to_add <- genes_present[!genes_present %in% data_lst[[f]]$V1]

  col_names <- colnames(data_lst[[f]])

  additional_rows <- data.frame(matrix(0,
    nrow = length(genes_to_add),
    ncol = length(col_names)
  ))
  colnames(additional_rows) <- col_names
  additional_rows$V1 <- genes_to_add

  # add the additional, empty probes and arrange in a common order
  dt_out <- data_lst[[f]] %>%
    bind_rows(additional_rows) %>%
    .[match(genes_present, .$V1), ]
  
  if (any(dt_out$V1 != genes_present)) {
    stop("Check order is correct")
  }

  
  
  file_write <- strsplit(f, "_?_(.*?)_?")[[1]][2] 
  empty_probes <- ! dt_out$V1 %in% genes_to_add
  empty_probes_dt[[file_write]] <- empty_probes
  
  file_write %<>% paste0(., ".csv")

  fwrite(dt_out, file = file_write, row.names = F)
}

fwrite(empty_probes_dt, file = "probes_present_per_dataset.csv")

summary(dt_out[, 1:5])

genes_present_dt <- data.table(Probe_ID = genes_present)
fwrite(genes_present_dt, "probe_IDs_present.csv")
