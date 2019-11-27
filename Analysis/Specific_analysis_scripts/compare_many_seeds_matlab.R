

library(tibble)
library(pheatmap)
library(purrr)
library(data.table)
library(png)
library(gridExtra)

# Function courtesy of https://sebastiansauer.github.io/two-plots-rmd/
# Combines 3 pngs (I extended it bravely)
comb3pngs <- function(imgs, bottom_text = NULL) {
  img1 <- grid::rasterGrob(as.raster(png::readPNG(imgs[1])),
                           interpolate = FALSE
  )
  img2 <- grid::rasterGrob(as.raster(png::readPNG(imgs[2])),
                           interpolate = FALSE
  )
  img3 <- grid::rasterGrob(as.raster(png::readPNG(imgs[3])),
                           interpolate = FALSE
  )
  gridExtra::grid.arrange(img1, img2, img3, ncol = 3, bottom = bottom_text)
}

many_seeds_5000_tbl_filename <- paste0(
  "/home/MINTS/sdc56/Desktop/MDI_small_geneset_outputs/",
  "Many_seeds_small_7_output_5000/compare_tibble.rds"
)

matlab_no_norm_445_tbl_filename <- paste0(
  "/home/MINTS/sdc56/Desktop/MDI_small_geneset_outputs/",
  "matlab_output_no_vsn_no_norm/compare_tibble.rds"
)

probe_key_file <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/probe_key.csv"

plot_type <- ".png"

many_seeds_5000_tbl <- readRDS(many_seeds_5000_tbl_filename)
matlab_no_norm_445_tbl <- readRDS(matlab_no_norm_445_tbl_filename)

probe_key <- data.table::fread(probe_key_file)

probe_names <- colnames(many_seeds_5000_tbl$mdi_allocation[1][[1]])

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

save_dir <- "/home/MINTS/sdc56/Desktop/MDI_small_geneset_outputs/Comparison/"
ph_names <- c("many_seeds", "matlab", "expression_data")


for (dataset in datasets) {
  # dataset <- "CD14"
  # Directory for specific cell type
  loc_dir <- paste0(save_dir, dataset)
  dir.create(loc_dir, showWarnings = FALSE)

  # File names for the three heatmaps
  ph_filenames <- paste0(loc_dir, "/", ph_names, plot_type)

  # Pull out the relevant data from the tibble
  many_seeds_sim <- many_seeds_5000_tbl$similarity_matrix[many_seeds_5000_tbl$dataset == dataset][[1]] %>%
    magrittr::set_colnames(gene_id) %>%
    magrittr::set_rownames(gene_id)

  matlab_sim <- matlab_no_norm_445_tbl$similarity_matrix[matlab_no_norm_445_tbl$dataset == dataset][[1]] %>%
    magrittr::set_colnames(gene_id) %>%
    magrittr::set_rownames(gene_id)

  expression_data <- matlab_no_norm_445_tbl$expression_data[matlab_no_norm_445_tbl$dataset == dataset][[1]] %>%
    magrittr::set_rownames(gene_id)

  # Create the heatmap for the PSM for MDI based on multiple seeds
  ph_many_seeds <- pheatmap(many_seeds_sim,
    main = paste0(dataset, ": PSM for many seeds"),
    filename = ph_filenames[1],
    cellheight = 3.0
  )

  # Extract the ordering from this
  row_order <- ph_many_seeds$tree_row$order

  # Heatmap of PSM from MATLAB output
  ph_matlab <- pheatmap(matlab_sim[row_order, row_order],
    cluster_rows = F,
    cluster_cols = F,
    main = paste0(dataset, ": PSM for MATLAB"),
    filename = ph_filenames[2],
    cellheight = 3.0
  )

  # Heatmap of expression data
  ph_expresssion <- pheatmap(expression_data[row_order, ],
    cluster_rows = F,
    cluster_cols = F,
    main = paste0(dataset, ": expression data"),
    filename = ph_filenames[3],
    cellheight = 3.0
  )
  
  
  ph_png_filenames <- paste0(loc_dir, "/", ph_names, ".png")
  comp_save_title <- paste0(loc_dir, "/Comparison.pdf")
  
  # if (plot_type == ".pdf") {
    pdf(comp_save_title)
  # } else {
    # png(comp_save_title)
  # }
  comb3pngs(ph_png_filenames)
  dev.off()
}

# 
# 
# 
# png_home <- "/home/MINTS/sdc56/Desktop/MDI_small_geneset_outputs/Comparison/CD4/"
# png_files <- c("expression_data.png", "many_seeds.png", "matlab.png")
# 
# png_files <- paste0(png_home, png_files)
# 
# comb3pngs(png_files[c(2, 3, 1)])
# 
# # Define "fused" genes
# fusion_threshold <- 0.5
# count <- 0
# 
# num_combinations <- choose(num_datasets, 2)
# 
# fusion_prob_lst <- list()
# fused_probes <- list()
# fusion_tbl <- tibble(
#   probabilities = vector("list", num_combinations),
#   probes_ids = vector("list", num_combinations)
# )
# 
# matlab_non_zero_probes_1 <- rowSums(matlab_no_norm_445_tbl$expression_data[1][[1]]) != 0
# matlab_non_zero_probes_2 <- rowSums(matlab_no_norm_445_tbl$expression_data[2][[1]]) != 0
# 
# many_seeds_non_zero_probes_1 <- rowSums(many_seeds_5000_tbl$expression_data[1][[1]]) != 0
# many_seeds_non_zero_probes_2 <- rowSums(many_seeds_5000_tbl$expression_data[2][[1]]) != 0
# 
# matlab_fusion_tbl <- many_seeds_fusion_tbl <- fusion_tbl
# 
# for(i in 1:(num_datasets - 1)){
#   for(j in (i + 1) : num_datasets){
#     
#     count <- count + 1
#     
#     .fusion_count <- matlab_no_norm_445_tbl$mdi_allocation[i][[1]] == matlab_no_norm_445_tbl$mdi_allocation[j][[1]]
#     
#     matlab_fusion_tbl$probabilities[count][[1]] <- .fusion_prob <- (1/nrow(.fusion_count)) * colSums(.fusion_count)
#     
#     matlab_fusion_tbl$probes_ids[count][[1]] <- .fused_probes_ind <- .fusion_prob > fusion_threshold
#     
#     .fusion_count <- many_seeds_5000_tbl$mdi_allocation[i][[1]] == many_seeds_5000_tbl$mdi_allocation[j][[1]]
#     
#     many_seeds_fusion_tbl$probabilities[count][[1]] <- .fusion_prob <- (1/nrow(fusion_count)) * colSums(fusion_count)
#     
#     many_seeds_fusion_tbl$probes_ids[count][[1]] <- .fused_probes_ind <- .fusion_prob > fusion_threshold
#   }
# }
# # 
# # for(i in 1:(num_datasets - 1)){
# #   for(j in (i + 1) : num_datasets){
# #     
# #     count <- count + 1
# #     
# #     .fusion_count <- matlab_no_norm_445_tbl$mdi_allocation[i][[1]] == matlab_no_norm_445_tbl$mdi_allocation[j][[1]]
# #     
# #     fusion_tbl$probabilities[count][[1]] <- .fusion_prob <- (1/nrow(fusion_count)) * colSums(fusion_count)
# #     
# #     fusion_tbl$probes_ids[count][[1]] <- .fused_probes_ind <- fusion_prob > fusion_threshold
# #   }
# # }
# 
# 
# ph1 <- pheatmap(many_seeds_5000_tbl$similarity_matrix[1][[1]])
# 
# pheatmap(matlab_no_norm_445_tbl$similarity_matrix[1][[1]][ph1$tree_row$order, ph1$tree_row$order], cluster_rows = 0, cluster_cols = 0)
# 
# ph_1 <- pheatmap(matlab_no_norm_445_tbl$similarity_matrix[1][[1]][fused_probes_ind, fused_probes_ind])
# # 
# pheatmap(matlab_no_norm_445_tbl$similarity_matrix[2][[1]][fused_probes_ind, fused_probes_ind])
# 
# many_seeds_sim_1 <- many_seeds_5000_tbl$similarity_matrix[1][[1]]
# 
# many_seeds_sim_1[many_seeds_fused_probes_1, many_seeds_fused_probes_1] %>% 
#   pheatmap()
# 
# many_seeds_sim_1[many_seeds_non_zero_probes_1, many_seeds_non_zero_probes_1] %>% 
#   pheatmap()
# 
# many_seeds_fused_probes_1 <- many_seeds_fusion_tbl$probes_ids[1][[1]]
# fused_non_zero_probes <- many_seeds_fused_probes_1 & many_seeds_non_zero_probes_1
#   
# 
# many_seeds_sim_1[fused_non_zero_probes, fused_non_zero_probes] %>% 
#   pheatmap()
# 
# pheatmap(many_seeds_5000_tbl$similarity_matrix[1][[1]][many_seeds_fusion_tbl$probes_ids[1][[1]], many_seeds_fusion_tbl$probes_ids[1][[1]]])
# 
# pheatmap(many_seeds_5000_tbl$similarity_matrix[2][[1]][many_seeds_fusion_tbl$probes_ids[1][[1]], many_seeds_fusion_tbl$probes_ids[1][[1]]])
# 
# pheatmap(many_seeds_5000_tbl$similarity_matrix[1][[1]][many_seeds_5000_tbl$fused_non_zero_probes[1][[1]], many_seeds_5000_tbl$fused_non_zero_probes[1][[1]]])
# 
# many_seeds_5000_tbl$non_zero_probes
# many_seeds_5000_tbl$fused_probes
# 
# dataset_i <- "CD14"
# dataset_j <- "CD19"
# i <- j <- k <- 1
# start_index <- 1
# eff_n_iter <- 5000
# colSums(many_seeds_5000_tbl$mdi_allocation[many_seeds_5000_tbl$dataset == dataset_i][[k]] == many_seeds_5000_tbl$mdi_allocation[many_seeds_5000_tbl$dataset == dataset_j][[k]])
# 
# head(many_seeds_5000_tbl$mdi_allocation[1][[1]][1:4, 1:4])
# head(many_seeds_5000_tbl$mdi_allocation[2][[1]][1:4, 1:4])
# 
# %>%
#   magrittr::extract(start_index:eff_n_iter,)
# 
# ind_used <- matlab_no_norm_445_tbl$fused_probes[[1]]
# 
# 
# sum(matlab_fusion_tbl$probes_ids[[6]])
# sum(matlab_no_norm_445_tbl$fused_probes[[1]])
# ph_1 <- pheatmap(matlab_no_norm_445_tbl$similarity_matrix[1][[1]][matlab_fusion_tbl$probes_ids[[6]], matlab_fusion_tbl$probes_ids[[6]]])
# pheatmap(matlab_no_norm_445_tbl$similarity_matrix[[1]][ind_used, ind_used])
# 
# pheatmap(matlab_no_norm_445_tbl$similarity_matrix[1][[1]][matlab_fusion_tbl$probes_ids[[2]], matlab_fusion_tbl$probes_ids[[2]]])
# pheatmap(matlab_no_norm_445_tbl_2$similarity_matrix[[1]][matlab_no_norm_445_tbl_2$fused_non_zero_probes[[1]][[3]], matlab_no_norm_445_tbl_2$fused_non_zero_probes[[1]][[3]]])
