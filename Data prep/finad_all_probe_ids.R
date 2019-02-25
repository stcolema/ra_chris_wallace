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

# Read in data
# files_present <- list.files(path = args$dir)
# file_name <- grep(args$extension, files_present, value = TRUE) %>%
#   paste0(args$dir, "/", .)
setwd("~/Desktop/MDI/Data/My_end")
# setwd("~/Desktop/Na_filled_data/")
files_present <- list.files(path = "~/Desktop/MDI/Data/My_end")
# files_present <- list.files(path = "~/Desktop/Na_filled_data/")
file_name <- grep(".csv", files_present, value = TRUE)

do_pca <- F
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
  if (do_pca) {
    pca_lst[[f]] <- prcomp(data_lst[[f]][, -1])

    pca_plot_lst[[f]] <- fviz_pca_ind(pca_lst[[f]],
      geom.ind = "point", # show points only (nbut not "text")
      col.ind = "contrib"
    ) +
      scale_color_gradient2(
        low = "black", mid = "blue",
        high = "red", midpoint = 4
      )
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
if (do_pca) {
  print(pca_plot_lst$transposed_CD14_GE_Corrected4_Covars.csv)
  print(pca_plot_lst$transposed_IL_GE_Corrected4_Covars.csv)
  print(pca_plot_lst$transposed_CD19_GE_Corrected4_Covars.csv)
  print(pca_plot_lst$transposed_CD15_GE_Corrected4_Covars.csv)
  print(pca_plot_lst$transposed_CD4_GE_Corrected4_Covars.csv)
  print(pca_plot_lst$transposed_CD8_GE_Corrected4_Covars.csv)
  print(pca_plot_lst$transposed_PLA_GE_Corrected4_Covars.csv)
  print(pca_plot_lst$transposed_RE_GE_Corrected4_Covars.csv)
  print(pca_plot_lst$transposed_TR_GE_Corrected4_Covars.csv)
}
# length(genes_present)

# Move the genes present to a list
genes_present %<>% unname() %>% unlist()

empty_probes_dt <- data.table(matrix(NA,
                  nrow = length(genes_present),
                  ncol = length(file_name) + 1
))

files_to_write <- strsplit(names(data_lst), "_?_(.*?)_?") %>% 
  lapply("[[", 2) %>% 
  unlist()

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
