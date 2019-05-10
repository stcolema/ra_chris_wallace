

library(data.table)
library(pheatmap)

data_dir <- "/home/MINTS/sdc56/Desktop/Matlab_input_small_names"
gen_ph_title <- ": heatmap of expression data"
gen_ph_file_name <- "/home/MINTS/sdc56/Desktop/pheatmap_"
file_type <- ".png"
n_genes <- 102

mdi_input_files <- list.files(path = data_dir, full.names = T, include.dirs = F) %>%
  grep("csv", ., value = TRUE) %>% 
  .[-c(2, 7)]
num_files <- length(mdi_input_files)
file_names <- basename(tools::file_path_sans_ext(mdi_input_files))

datasets <- file_names %>% stringr::str_replace("_sma_mat_nv", "")

mega_df <- data.frame(matrix(nrow = n_genes, ncol = 0)) 



data_files <- list()
for(i in 1 : num_files){
  
  file_name <- gen_ph_file_name %>% 
    paste0(datasets[[i]], file_type)
  
  ph_title <- datasets[[i]] %>% 
    paste0(gen_ph_title)
  
  f <- mdi_input_files[[i]]
  data_files[[i]] <- fread(f)
  data_files[[i]][is.na(data_files[[i]])] <- 0
  data_files[[i]][, -1] %>% 
    pheatmap(filename = file_name, main = ph_title)
  
  mega_df <- mega_df %>% 
    dplyr::bind_cols(data_files[[i]][, -1])
  
}

row.names(mega_df) <- data_files[[1]][,1] %>% unlist()

big_file_name <- gen_ph_file_name %>% 
  paste0("All_datasets", file_type)

big_ph_title <- "All datasets" %>% 
  paste0(gen_ph_title)

pheatmap(mega_df, filename = big_file_name, main = big_ph_title)

