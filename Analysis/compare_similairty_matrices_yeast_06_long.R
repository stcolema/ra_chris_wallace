

function_dir <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/Analysis_script_functions/"

function_scripts <- c(
  "plot_similarity_matrices.R",
  "plot_comparison_expression_clustering.R",
  "create_psm.R"
)

for (f in paste0(function_dir, function_scripts)) {
  source(f)
}

# For tibbles
library(tibble) # for dataframe of lists

# For data wrangling
library(dplyr)

# Heatmapping
library(pheatmap) # install.packages("pheatmap", dep = T)

# Colour palettes
library(RColorBrewer)

# Find the MDI output directories
yeast_dir <- "~/Desktop/Yeast_output_3_datasets_long_diffuse"
yeast_output_dirs <- list.dirs(yeast_dir, recursive = F)

ind_to_drop_for_name <- nchar("/home/MINTS/sdc56/Desktop/Yeast_output_3_datasets_long_diffuse/Diffuse_long_06_")

yeast_comparison_dir <- "~/Desktop/Comparison_yeast_3_datasets/"
dir.create(yeast_comparison_dir, showWarnings = FALSE)

tibble_name <- "compare_tibble.rds" %>% 
  paste0(yeast_output_dirs, "/", .)

n_files <- length(yeast_output_dirs)

# imnputs for heatmaps
col_pal_sim <- colorRampPalette(c("white", "#146EB4"))(100)
palette_length <- length(col_pal_sim)

sim_breaks <- c(
  seq(1 / palette_length, 1, length.out = palette_length)
)

col_pal_expr <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(100)

expr_breaks <- define_breaks(col_pal_expr)

# Read in the tibbles
tib_list <- lapply(tibble_name, readRDS)

# The datasets present
datasets <- tib_list[[1]]$dataset

dataset_dirs <- paste0(yeast_comparison_dir, "/", datasets)
lapply(dataset_dirs, dir.create, showWarnings = F)

# Vector of relevant strings for differentiating save files
save_strings <- yeast_output_dirs %>% 
  sapply(function(x){
    substr(x, ind_to_drop_for_name + 1, nchar(x))
  }) %>% unname()

# Similar to above but nicer for titles
title_strings <- save_strings %>% 
  str_replace_all("_"," ") %>%
  str_to_sentence


for(d in datasets){

  curr_ind <- which(tib_list[[1]]$dataset == d)
  curr_save_dir <- dataset_dirs[curr_ind]
  
for(i in 1:(n_files - 2)){
  
  save_str_i <- save_strings[i]
  title_str_i <- title_strings[i]
  
  sim_i <- tib_list[[i]]$similarity_matrix[[curr_ind]]
  
  row_order <- hclust(dist(sim_i))$order
  
  sim_i <- sim_i[row_order, row_order]
  
  ph_sim_i <- pheatmap(sim_i,
     cluster_rows = F,
     cluster_cols = F,
     color = col_pal_sim,
     breaks = sim_breaks,
     silent = T
  )$gtable
  
  for(j in (i+1):(n_files - 1)){
    
    save_str_j <- save_strings[j]
    title_str_j <- title_strings[j]
    
    sim_j <- tib_list[[j]]$similarity_matrix[[curr_ind]][row_order, row_order]
    
    ph_sim_j <- pheatmap(sim_j,
                       cluster_rows = F,
                       cluster_cols = F,
                       color = col_pal_sim,
                       breaks = sim_breaks,
                       silent = T
    )$gtable
    
    for(k in (j+1):n_files){
      
      save_str_k <- save_strings[k]
      title_str_k <- title_strings[k]
      
      sim_k <- tib_list[[k]]$similarity_matrix[[curr_ind]][row_order, row_order]
      
      ph_sim_k <- pheatmap(sim_k,
                         cluster_rows = F,
                         cluster_cols = F,
                         color = col_pal_sim,
                         breaks = sim_breaks,
                         silent = T
      )$gtable
      
      loc_save_name <- paste(save_str_i, save_str_j, save_str_k, sep = "_")
      save_name <- paste0(curr_save_dir, "/", loc_save_name, ".png")

      curr_title <- paste0(d,
                           ": Similarity matrices for ",
                           title_str_i,
                           " (defines order), ",
                           title_str_j,
                           ", and ",
                           title_str_k)
      
      combine_pheatmaps(list(ph_sim_i, ph_sim_j, ph_sim_k),
                        save_name = save_name,
                        main = curr_title
      )
      
    }
  }
}

}
