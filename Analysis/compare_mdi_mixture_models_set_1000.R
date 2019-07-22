

library(data.table)
library(pheatmap)
library(cowplot)
library(ggplot2)
library(magrittr)

source("~/Desktop/ra_chris_wallace/Analysis/Analysis_script_functions/plot_comparison_expression_clustering.R")

tib_all <- readRDS("~/Desktop/Final_set//CD/Analysis/compare_tibble.rds")


tissue_type <- c("CD", "Colon")


for(tissue in tissue_type){
  all_dir <- paste0("~/Desktop/Final_set/", tissue, "/Analysis/")
  specific_dir <- "~/Desktop/Final_set/Dataset_specific/"
  
  gen_tib_name <- "compare_tibble.rds"
  save_dir <- "~/Desktop/Final_set/Comparison_mdi_mixture_model/"
  gen_ph_name <- paste0("comparison_", tissue, "_specific_sim_cor.png")

  

tib_all <- readRDS(paste0(all_dir, gen_tib_name)) 
datasets <- tib_all$dataset
n_datasets <- length(datasets)



col_pal <- colorRampPalette(c("#FF9900", "white", "#146EB4"))(100)
my_breaks <- define_breaks(col_pal)

cor_pal <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(100)
cor_breaks <- define_breaks(cor_pal)

# tib_specific <- vector("list", n_datasets) %>% set_names(datasets)


for(d in datasets){
  
  if(d == "CD15") {
    next
  }
  
  curr_ind <- which(tib_all$dataset == d)
  # tib_specific[[d]] <- 
  tib_specific <-  readRDS(paste0(specific_dir, d, "/", gen_tib_name))
  
  all_sim <- tib_all$similarity_matrix[[curr_ind]]
  spec_sim <- tib_specific$similarity_matrix[[1]]
  cor_mat <- tib_all$correlation_matrix[[curr_ind]]
  
  row_order <- hclust(dist(all_sim))$order
  
  # Re order the matrices to have a common row order
  all_sim <- all_sim[row_order, row_order]
  spec_sim <- spec_sim[row_order, row_order]
  cor_mat <- cor_mat[row_order, row_order]
  
  
  ph_list <- list()
  
  ph_title <- paste0(d, ": Comparison MDI to mixure model")
  ph_save <- paste0(save_dir, d, "_", gen_ph_name)
  
  # Create a heatmap of the data without clustering
  ph_list[[1]] <- pheatmap(all_sim,
                           cluster_rows = F,
                           cluster_cols = F,
                           color = col_pal,
                           breaks = my_breaks,
                           show_rownames = F,
                           show_colnames = F,
                           silent = TRUE
  )$gtable
  
  ph_list[[2]] <- pheatmap(spec_sim,
                           cluster_rows = F,
                           cluster_cols = F,
                           color = col_pal,
                           breaks = my_breaks,
                           show_rownames = F,
                           show_colnames = F,
                           silent = TRUE
  )$gtable
  
  ph_list[[3]] <- pheatmap(cor_mat,
                           cluster_rows = F,
                           cluster_cols = F,
                           color = cor_pal,
                           breaks = cor_breaks,
                           show_rownames = F,
                           show_colnames = F,
                           silent = TRUE
  )$gtable
  
  # Combine these in a grid and save
  combine_pheatmaps(ph_list, save_name = ph_save, main = ph_title, font_size = 20)
  
}
}
