


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

find_file_names <- function(weights, 
                            n_weights, 
                            iter,
                            long_runs_dir_parts,
                            many_seeds_500_dir_parts,
                            gen_tibble_name = "compare_tibble.rds"){
  
  files_of_interest <- vector("list", n_weights)
  names(files_of_interest) <- paste0("Weights_", weights)
  
  for(i in 1 : n_weights){
    long_runs_dir <- paste0(long_runs_dir_parts[1], weights[i], long_runs_dir_parts[2], iter, "/")
    many_seeds_500_dir <- paste0(many_seeds_500_dir_parts, weights[i], "/")
    all_dirs <- c(many_seeds_500_dir, long_runs_dir)
    
    gen_tibble_name <- "compare_tibble.rds"
    
    all_files <- paste0(yeast_dir, all_dirs, gen_tibble_name)
    
    files_of_interest[[i]] <- all_files
  }
  files_of_interest
}


col_pal <- colorRampPalette(c("white", "#146EB4"))(100)
palette_length <- length(col_pal)

breaks <- c(
  seq(1 / palette_length, 1, length.out = palette_length)
)



yeast_dir <- "/home/MINTS/sdc56/Desktop/Yeast_MDI/Yeast_output/Diffuse/"

long_runs_dir_parts <- c("/Long_runs/iter_", "_seed_")

many_seeds_500_dir_parts <- "500_iter_"

n_versions_many_seeds <- 1

weights <- 1:5
n_weights <- length(weights)

iter <- 1:5

n_levels <- n_versions_many_seeds + (length(iter))

# long_runs_dir <- paste0(long_runs_dir_parts[1], weights[1], long_runs_dir_parts[2], iter, "/")
# many_seeds_500_dir <- paste0(many_seeds_500_dir_parts, weights, "/")
# all_dirs <- c(long_runs_dir, many_seeds_500_dir)
# 
# gen_tibble_name <- "compare_tibble.rds"
# 
# all_files <- paste0(yeast_dir, all_dirs, gen_tibble_name)

files_of_interest <- find_file_names(weights, 
  n_weights, 
  iter,
  long_runs_dir_parts,
  many_seeds_500_dir_parts,
  gen_tibble_name = "compare_tibble.rds"
  )

source_files <- c( "many_500", paste0("long_", iter))
my_tibble <- tibble(Weight = numeric(),
                        Source = character(),
                        Similarity_matrix = list(),
                        Dataset = character(),
                        Expression_data = list())

rel_for_ave <- paste0("long_", iter)

for(i in 1 : n_weights){
  
  w <- weights[i]
  
  for(j in 1 : n_levels){
    
    f <-  files_of_interest[[i]][[j]]
    
    s <- source_files[j]
    t <- readRDS(f)
    sim <- t$similarity_matrix
    d <- t$dataset
    e_d <- t$expression_data
    
    my_tib <- tibble(Weight = w, Source = s, Similarity_matrix = sim, Dataset = d, Expression_data = e_d)
    
    my_tibble <- bind_rows(my_tibble, my_tib)
    
  }
  
  for(j in 1:n_datasets){
    
    d <- datasets[j]
    
    ave_sim_ind <- which(my_tibble$Weight == w 
                         & my_tibble$Source %in% rel_for_ave
                         & my_tibble$Dataset == d)
    sim_list <- my_tibble$Similarity_matrix[ave_sim_ind]
    
    e_d <- my_tibble$Expression_data[ave_sim_ind][[1]]
    
    n_sims <- length(ave_sim_ind)
    
    ave_sim <- Reduce(`+`, sim_list) / n_sims
    
    my_tib <- tibble(Weight = w, 
                     Source = "average_long",
                     Similarity_matrix = list(ave_sim), 
                     Dataset = d, 
                     Expression_data = list(e_d)
                     )
    my_tibble <- bind_rows(my_tibble, my_tib)
  }
  
  
}
# 
# 
# # n_files <- length(all_files)
# # 
# # for(i in 1 : n_levels){
# #   f <-  all_files[[i]]
# #   my_tibble_list[[i]] <- readRDS(all_files[[i]])
# # }
# 
# file_names <- 
# 
# sim_list_j <- list()
# 
# datasets <-  unique(my_tibble_list[[1]]$dataset)
# n_datasets <- length(datasets)
# 
# sim_list_list <- list()
# 
# for(j in 1 : n_datasets) {
#   sim_list_list[[j]] <- list()
#   curr_dataset <- datasets[j]
#   curr_index <- which(my_tibble_list[[i]]$dataset == curr_dataset)
#   
#   for(i in 1 : n_levels){
#     
#     sim_list_list[[j]][[i]] <- my_tibble_list[[i]]$similarity_matrix[[curr_index]]
#     
#   }
#   
# }


source_files_update <- my_tibble$Source %>% unique()
n_levels_update <- n_levels + 1

gen_save_name <- "/home/MINTS/sdc56/Desktop/Yeast_MDI/Yeast_output/Diffuse/Comparison_diffuse/" %T>%
  dir.create(showWarnings = F)

gen_file_name <- "Sim_mat_comparison_"
file_type <- ".png"

print_comparison_plots(my_tibble, weights, datasets, source_files_update, n_levels_update,
  gen_save_name = "/home/MINTS/sdc56/Desktop/Yeast_MDI/Yeast_output/Diffuse/Comparison_diffuse/", 
  gen_file_name = "Sim_mat_comparison_",
  file_type = ".png",
  col_pal = colorRampPalette(c("white", "#146EB4"))(100),
  breaks = c(seq(1 / length(col_pal), 1, length.out = length(col_pal)))
)



# for(w in weights){
#   w_dir <- paste0(gen_save_name, w, "/") %T>%
#     dir.create(showWarnings = F)
#   for(d in datasets){
#     
#     d_dir <-  paste0(gen_save_name, w, "/", d, "/") %T>%
#       dir.create(showWarnings = F)
#     
#     for(i in 1:(n_levels_update - 2)){
#       d_i <- source_files_update[i]
#       ind_i <- which(my_tibble$Source == d_i 
#                      & my_tibble$Dataset == d
#                      & my_tibble$Weight == w)
#       
#       sim_i <- my_tibble$Similarity_matrix[[ind_i]]
#       
#       for(j in (i + 1):(n_levels_update - 1)){
#         d_j <- source_files_update[j]
#         ind_j <- which(my_tibble$Source == d_j
#                        & my_tibble$Dataset == d
#                        & my_tibble$Weight == w)
#         
#         sim_j <- my_tibble$Similarity_matrix[[ind_j]]
#         
#         for(k in (j + 1): n_levels_update){
#           d_k <- source_files_update[k]
#           ind_k <- which(my_tibble$Source == d_k
#                          & my_tibble$Dataset == d
#                          & my_tibble$Weight == w)
#           
#           sim_k <- my_tibble$Similarity_matrix[[ind_k]]
#           
#           ph_title <- paste0(d, 
#                              ": Comparison of clustering for ",
#                              d_i,
#                              ", ",
#                              d_j,
#                              " and ",
#                              d_k)
#           
#           save_name <- paste0(
#             d_dir,
#             gen_file_name,
#             d_i,
#             "_",
#             d_j,
#             "_",
#             d_k,
#             file_type
#           )
#           
#           heatmap_wrapper_sim_expr_corr(sim_i, sim_j, sim_k,
#                                         ph_title = ph_title,
#                                         save_name = save_name,
#                                         col_pal_sim = col_pal,
#                                         col_pal_expr = col_pal,
#                                         expr_breaks = breaks,
#                                         sim_breaks = breaks,
#                                         font_size = 20,
#                                         expr_col_order = F)
#           
#         }
#       }
#     }
#   }
# }

print_comparison_plots <- function(my_tibble, weights, datasets, source_files, n_levels,
                                   gen_save_name = "/home/MINTS/sdc56/Desktop/Yeast_MDI/Yeast_output/Diffuse/Comparison_diffuse/", 
                                   gen_file_name = "Sim_mat_comparison_",
                                   file_type = ".png",
                                   col_pal = colorRampPalette(c("white", "#146EB4"))(100),
                                   breaks = c(seq(1 / length(col_pal), 1, length.out = length(col_pal)))){

  dir.create(gen_save_name, showWarnings = F)
  
for(w in weights){
  w_dir <- paste0(gen_save_name, w, "/") %T>%
    dir.create(showWarnings = F)
  for(d in datasets){
    
    d_dir <-  paste0(gen_save_name, w, "/", d, "/") %T>%
      dir.create(showWarnings = F)
    
    for(i in 1:(n_levels - 2)){
      d_i <- source_files[i]
      ind_i <- which(my_tibble$Source == d_i 
                     & my_tibble$Dataset == d
                     & my_tibble$Weight == w)
      
      sim_i <- my_tibble$Similarity_matrix[[ind_i]]
      
      for(j in (i + 1):(n_levels - 1)){
        d_j <- source_files[j]
        ind_j <- which(my_tibble$Source == d_j
                       & my_tibble$Dataset == d
                       & my_tibble$Weight == w)
        
        sim_j <- my_tibble$Similarity_matrix[[ind_j]]
        
        for(k in (j + 1): n_levels){
          d_k <- source_files[k]
          ind_k <- which(my_tibble$Source == d_k
                         & my_tibble$Dataset == d
                         & my_tibble$Weight == w)
          
          sim_k <- my_tibble$Similarity_matrix[[ind_k]]
          
          ph_title <- paste0(d, 
                             ": Comparison of clustering for ",
                             d_i,
                             ", ",
                             d_j,
                             " and ",
                             d_k)
          
          save_name <- paste0(
            d_dir,
            gen_file_name,
            d_i,
            "_",
            d_j,
            "_",
            d_k,
            file_type
          )
          
          heatmap_wrapper_sim_expr_corr(sim_i, sim_j, sim_k,
                                        ph_title = ph_title,
                                        save_name = save_name,
                                        col_pal_sim = col_pal,
                                        col_pal_expr = col_pal,
                                        expr_breaks = breaks,
                                        sim_breaks = breaks,
                                        font_size = 20,
                                        expr_col_order = F)
          
        }
      }
    }
  }
}

}

# names(my_tibble_list) <- all_dirs
# 
# sim_500_1_1 <- my_tibble_list$`500_iter_1/`$similarity_matrix[[1]]
# sim_500_2_1 <- my_tibble_list$`500_iter_2/`$similarity_matrix[[1]]
# sim_500_3_1 <- my_tibble_list$`500_iter_3/`$similarity_matrix[[1]]
# sim_500_4_1 <- my_tibble_list$`500_iter_4/`$similarity_matrix[[1]]
# sim_500_5_1 <- my_tibble_list$`500_iter_5/`$similarity_matrix[[1]]
# sim_long_1_1_1 <- my_tibble_list$`/Long_runs/iter_1_seed_1/`$similarity_matrix[[1]]
# sim_long_1_2_1 <- my_tibble_list$`/Long_runs/iter_1_seed_2/`$similarity_matrix[[1]]
# sim_long_1_3_1 <- my_tibble_list$`/Long_runs/iter_1_seed_3/`$similarity_matrix[[1]]
# sim_long_1_4_1 <- my_tibble_list$`/Long_runs/iter_1_seed_4/`$similarity_matrix[[1]]
# sim_long_1_5_1 <- my_tibble_list$`/Long_runs/iter_1_seed_5/`$similarity_matrix[[1]]
# avg_long_1_1 <- (sim_long_1_1_1 
#                  + sim_long_1_2_1
#                  + sim_long_1_3_1
#                  + sim_long_1_4_1 
#                  + sim_long_1_5_1) /5
# 
# sim_500_1_2 <- my_tibble_list$`500_iter_1/`$similarity_matrix[[5]]
# sim_500_2_2 <- my_tibble_list$`500_iter_2/`$similarity_matrix[[5]]
# sim_500_3_2 <- my_tibble_list$`500_iter_3/`$similarity_matrix[[5]]
# sim_500_4_2 <- my_tibble_list$`500_iter_4/`$similarity_matrix[[5]]
# sim_500_5_2 <- my_tibble_list$`500_iter_5/`$similarity_matrix[[5]]
# sim_long_1_1_2 <- my_tibble_list$`/Long_runs/iter_1_seed_1/`$similarity_matrix[[5]]
# sim_long_1_2_2 <- my_tibble_list$`/Long_runs/iter_1_seed_2/`$similarity_matrix[[5]]
# sim_long_1_3_2 <- my_tibble_list$`/Long_runs/iter_1_seed_3/`$similarity_matrix[[5]]
# sim_long_1_4_2 <- my_tibble_list$`/Long_runs/iter_1_seed_4/`$similarity_matrix[[5]]
# sim_long_1_5_2 <- my_tibble_list$`/Long_runs/iter_1_seed_5/`$similarity_matrix[[5]]
# avg_long_1_2 <- (sim_long_1_1_2
#                  + sim_long_1_2_2
#                  + sim_long_1_3_2
#                  + sim_long_1_4_2 
#                  + sim_long_1_5_2) /5
# 
# expr_data_1 <- my_tibble_list$`500_iter_1/`$expression_data[[1]]
# expr_data_2 <- my_tibble_list$`500_iter_1/`$expression_data[[5]]
# 
# heatmap_wrapper_sim_expr_corr(sim_500_1_1, expr_data_1, sim_long_1_1_1,
#                               ph_title = "Comparison consesnsus to long for 2 seeds",
#                               save_name = "/home/MINTS/sdc56/Desktop/yeast_diffuse_1_1_1.png",
#                               col_pal_sim = col_pal,
#                               col_pal_expr = col_pal,
#                               expr_breaks = breaks,
#                               sim_breaks = breaks,
#                               font_size = 13,
#                               expr_col_order = T)
# 
# heatmap_wrapper_sim_expr_corr(sim_500_1_1, expr_data_1, sim_long_1_2_1,
#                               ph_title = "Comparison consesnsus to long for 2 seeds",
#                               save_name = "/home/MINTS/sdc56/Desktop/yeast_diffuse_1_2_1.png",
#                               col_pal_sim = col_pal,
#                               col_pal_expr = col_pal,
#                               expr_breaks = breaks,
#                               sim_breaks = breaks,
#                               font_size = 13,
#                               expr_col_order = T)
# 
# heatmap_wrapper_sim_expr_corr(sim_500_1_1, expr_data_1, sim_long_1_3_1,
#                               ph_title = "Comparison consesnsus to long for 2 seeds",
#                               save_name = "/home/MINTS/sdc56/Desktop/yeast_diffuse_1_3_1.png",
#                               col_pal_sim = col_pal,
#                               col_pal_expr = col_pal,
#                               expr_breaks = breaks,
#                               sim_breaks = breaks,
#                               font_size = 13,
#                               expr_col_order = T)
# 
# heatmap_wrapper_sim_expr_corr(sim_500_1_1, expr_data_1, sim_long_1_4_1,
#                               ph_title = "Comparison consesnsus to long for 2 seeds",
#                               save_name = "/home/MINTS/sdc56/Desktop/yeast_diffuse_1_4_1.png",
#                               col_pal_sim = col_pal,
#                               col_pal_expr = col_pal,
#                               expr_breaks = breaks,
#                               sim_breaks = breaks,
#                               font_size = 13,
#                               expr_col_order = T)
# 
# heatmap_wrapper_sim_expr_corr(sim_500_1_1, expr_data_1, sim_long_1_5_1,
#                               ph_title = "Comparison consesnsus to long for 2 seeds",
#                               save_name = "/home/MINTS/sdc56/Desktop/yeast_diffuse_1_5_1.png",
#                               col_pal_sim = col_pal,
#                               col_pal_expr = col_pal,
#                               expr_breaks = breaks,
#                               sim_breaks = breaks,
#                               font_size = 13,
#                               expr_col_order = T)
# 
# heatmap_wrapper_sim_expr_corr(sim_500_1_1, expr_data_1, avg_long_1_1,
#                               ph_title = "Comparison consesnsus to long for 2 seeds",
#                               save_name = "/home/MINTS/sdc56/Desktop/yeast_diffuse_1_1_avg.png",
#                               col_pal_sim = col_pal,
#                               col_pal_expr = col_pal,
#                               expr_breaks = breaks,
#                               sim_breaks = breaks,
#                               font_size = 13,
#                               expr_col_order = T)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# heatmap_wrapper_sim_expr_corr(sim_500_1_2, expr_data_2, sim_long_1_1_2,
#                               ph_title = "Comparison consesnsus to long for 2 seeds",
#                               save_name = "/home/MINTS/sdc56/Desktop/yeast_diffuse_1_1_2.png",
#                               col_pal_sim = col_pal,
#                               col_pal_expr = col_pal,
#                               expr_breaks = breaks,
#                               sim_breaks = breaks,
#                               font_size = 13,
#                               expr_col_order = T)
# 
# heatmap_wrapper_sim_expr_corr(sim_500_1_2, expr_data_2, sim_long_1_2_2,
#                               ph_title = "Comparison consesnsus to long for 2 seeds",
#                               save_name = "/home/MINTS/sdc56/Desktop/yeast_diffuse_1_2_2.png",
#                               col_pal_sim = col_pal,
#                               col_pal_expr = col_pal,
#                               expr_breaks = breaks,
#                               sim_breaks = breaks,
#                               font_size = 13,
#                               expr_col_order = T)
# 
# heatmap_wrapper_sim_expr_corr(sim_500_1_2, expr_data_2, sim_long_1_3_2,
#                               ph_title = "Comparison consesnsus to long for 2 seeds",
#                               save_name = "/home/MINTS/sdc56/Desktop/yeast_diffuse_1_3_2.png",
#                               col_pal_sim = col_pal,
#                               col_pal_expr = col_pal,
#                               expr_breaks = breaks,
#                               sim_breaks = breaks,
#                               font_size = 13,
#                               expr_col_order = T)
# 
# heatmap_wrapper_sim_expr_corr(sim_500_1_2, expr_data_2, sim_long_1_4_2,
#                               ph_title = "Comparison consesnsus to long for 2 seeds",
#                               save_name = "/home/MINTS/sdc56/Desktop/yeast_diffuse_1_4_2.png",
#                               col_pal_sim = col_pal,
#                               col_pal_expr = col_pal,
#                               expr_breaks = breaks,
#                               sim_breaks = breaks,
#                               font_size = 13,
#                               expr_col_order = T)
# 
# heatmap_wrapper_sim_expr_corr(sim_500_1_2, expr_data_2, sim_long_1_5_2,
#                               ph_title = "Comparison consesnsus to long for 2 seeds",
#                               save_name = "/home/MINTS/sdc56/Desktop/yeast_diffuse_1_5_2.png",
#                               col_pal_sim = col_pal,
#                               col_pal_expr = col_pal,
#                               expr_breaks = breaks,
#                               sim_breaks = breaks,
#                               font_size = 13,
#                               expr_col_order = T)
# 
# heatmap_wrapper_sim_expr_corr(sim_500_1_2, expr_data_2, avg_long_1_2,
#                               ph_title = "Comparison consesnsus to long for 2 seeds",
#                               save_name = "/home/MINTS/sdc56/Desktop/yeast_diffuse_1_ave_2.png",
#                               col_pal_sim = col_pal,
#                               col_pal_expr = col_pal,
#                               expr_breaks = breaks,
#                               sim_breaks = breaks,
#                               font_size = 13,
#                               expr_col_order = T)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 

# heatmap_wrapper_sim_expr_corr(sim_500_1, sim_500_4, sim_500_5,
#                               ph_title = "Datasets 1, 4 and 5 for 1,000 seeds for 500 iterations",
#                               # save_name = "/home/MINTS/sdc56/Desktop/yeast_500_1_4_5.png",
#                               col_pal_sim = col_pal,
#                               col_pal_expr = col_pal,
#                               expr_breaks = breaks,
#                               sim_breaks = breaks,
#                               font_size = 13,
#                               expr_col_order = F)

