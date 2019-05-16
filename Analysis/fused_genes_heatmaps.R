


library(tidyverse)
library(magrittr)
library(pheatmap)

file_path <- "/home/MINTS/sdc56/Desktop/MDI_small_geneset_outputs/matlab_output_no_vsn_no_norm/"

compare_tibble_filename <- paste0(file_path, "compare_tibble.rds")
compare_tibble <- readRDS(compare_tibble_filename)

probe_key_file <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/probe_key.csv"

plot_type <- ".png"

probe_key <- data.table::fread(probe_key_file)

probe_names <- colnames(compare_tibble$mdi_allocation[1][[1]])

probe_key_rel <- probe_key[probe_key$ProbeID %in% probe_names, ]

# Pull out the Gene IDs in the correct order
gene_id <- probe_key_rel %>%
  .[match(probe_names, .$ProbeID)] %>%
  .$Gene

datasets <- c(
  "CD14",
  "CD19",
  "CD4",
  "CD8",
  "IL",
  "RE",
  "TR"
)

num_datasets <- length(datasets)

count <- list()
for(i in 1 : num_datasets){
  count[[i]] <- compare_tibble$fused_probes[[i]] %>% 
    lapply(sum) %>% 
    unlist() %>% 
    sum() %>% 
    `-`(102)
}


dir_name <- paste0(file_path, "Fusion_expression_data/")
dir.create(dir_name, showWarnings = FALSE)
generic_ph_title <- paste0(dir_name, "heatmap_")

for(i in 1:(num_datasets- 1)){
  d_i <- datasets[[i]]
  expression_data_i <- compare_tibble$expression_data[[i]]
  for(j in (i+ 1):num_datasets){
    fused_ind <- compare_tibble$fused_probes[[i]][[j]] # fused_non_zero_probes[[1]][[2]]
    non_fused_ind <- ! fused_ind 
    d_j <- datasets[[j]]
    
    expression_data_j <- compare_tibble$expression_data[[j]]
    
    ph1 <- pheatmap(expression_data_i[fused_ind, ])
    ph2 <- pheatmap(expression_data_j[fused_ind, ])
    
    col_order_1 <- ph1$tree_col$order
    col_order_2 <- ph2$tree_col$order
    
    new_expression_data <- cbind(expression_data_i[, col_order_1], expression_data_j[, col_order_2]) %>% 
      magrittr::set_rownames(gene_id)
    
    fused_ph_file_name <- paste0(generic_ph_title,
                                 "fused_genes_", 
                                 d_i,
                                 "_",
                                 d_j, 
                                 plot_type
                                 )
    
    pheatmap(new_expression_data[fused_ind, ],
             cluster_cols = F,
             gaps_col = length(col_order_1),
             main = paste("Fused probes for", d_i, "and", d_j),
             filename = fused_ph_file_name
    )
    
    unfused_ph_file_name <- paste0(generic_ph_title,
                                 "unfused_genes_", 
                                 d_i,
                                 "_",
                                 d_j, 
                                 plot_type
    )
    
    pheatmap(new_expression_data[non_fused_non_zeros, ],
             cluster_cols = F, 
             gaps_col = length(col_order_1),
             main = paste("Unfused probes for", d_i, "and", d_j),
             filename = unfused_ph_file_name
    )
    
  }
}