


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

gen_dir <- "/home/MINTS/sdc56/Desktop/Gen_data_output/"
gen_long_runs_dir <- c("Long", "Longer")
gen_many_seeds_dir <- "Many_seeds_"

tib_name <- "compare_tibble.rds"

types <- c("Longer", "500", "1000", "5000", "10000", "Long")
long_run_seeds <- c(9, 1:8, 10)

n_weights <- length(weights)
n_types <- length(types)

datasets <- c("MDItestdata1", "MDItestdata2", "MDItestdata3")

col_pal_sim <- colorRampPalette(c("white", "#146EB4"))(100)
palette_length <- length(col_pal_sim)

sim_breaks <- c(
  seq(1 / palette_length, 1, length.out = palette_length)
)


col_pal_corr <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(100)
corr_breaks <- define_breaks(col_pal_corr)

def_colours <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
col_pal_expr <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(100)
expr_breaks <- NA

sim_comparison_tibble <- tibble(
  Dataset = character(),
  Similarity_matrix = list(),
  Correlation_matrix = list(),
  Expression_data = list(),
  Similarity_matrix_ph = list(),
  Correlation_matrix_ph = list(),
  Expression_data_ph = list(),
  Sim_row_order = list(),
  Expr_col_order = list(),
  Type = character(),
  Seed = numeric()
)

n_genes <- 400

# Annotation labels based on generating method
true_clusters <- data.frame(
  Labels = as.factor(c(rep(1, 25), rep(2, 50), rep(3, 75), rep(4, 100), rep(5, 150)))
) %>% set_rownames(paste0("Gene_", 1:n_genes))

row_order <- order(true_clusters$Labels)

type_defining_order <- "NONE"
first_type <- which(types %in% type_defining_order)
types <- types[c(which(types %in% type_defining_order), which(! types %in% type_defining_order))]
seed_defing_order <- NA 


for (d in datasets) {
  # The relevant index within said tibble
  .curr_index <- which(datasets == d)
  
  .curr_expr <- read.csv(paste0(gen_dir, "/Data/", d, ".csv"), header = T, row.names = "V1") 
  
  # .curr_expr <- curr_tib$expression_data[[.curr_index]]
  
  expr_hc <- hclust(dist(t(.curr_expr)))
  expr_col_order <- expr_hc$order
  
  .curr_expr <- .curr_expr[row_order, expr_col_order]
  
  ph_expr <- pheatmap(.curr_expr,
                      cluster_rows = F,
                      cluster_cols = F,
                      color = col_pal_expr,
                      breaks = expr_breaks,
                      annotation_row = true_clusters,
                      labels_col = "",
                      silent = T
  )$gtable
  
  .curr_corr <- cor(t(.curr_expr))
  
  ph_corr <- pheatmap(.curr_corr,
                      cluster_rows = F,
                      cluster_cols = F,
                      color = col_pal_corr,
                      breaks = corr_breaks,
                      labels_col = "",
                      silent = T
  )$gtable
  
  for (j in 1:n_types) {
    curr_type <- types[j]
    
    
    if (curr_type == "Long" | curr_type == "Longer") {
      # for(curr_type in c("Long", "Longer")){
      for (k in long_run_seeds) {
        curr_dir <- paste0(gen_dir, curr_type, "_", k, "/")
        
        curr_tib <- readRDS(paste0(curr_dir, tib_name))
        
        
        .curr_index <- which(curr_tib$dataset == d)
        
        .curr_sim <- curr_tib$similarity_matrix[[.curr_index]] %>% 
          extract(row_order, row_order)
        
        .curr_corr <- curr_tib$correlation_matrix[[.curr_index]][row_order, row_order]
        .curr_expr <- curr_tib$expression_data[[.curr_index]][row_order, expr_col_order]
        
        ph_sim <- pheatmap(.curr_sim,
                           cluster_rows = F,
                           cluster_cols = F,
                           color = col_pal_sim,
                           breaks = sim_breaks,
                           annotation_row = true_clusters,
                           labels_col = "",
                           silent = T
        )$gtable
        
        # Add to a tibble
        sim_comparison_tibble_entry <- tibble(
          Dataset = d,
          Similarity_matrix = list(.curr_sim),
          Correlation_matrix = list(.curr_corr),
          Expression_data = list(.curr_expr),
          Similarity_matrix_ph = list(ph_sim),
          Correlation_matrix_ph = list(ph_corr),
          Expression_data_ph = list(ph_expr),
          Sim_row_order = list(row_order),
          Expr_col_order = list(expr_col_order),
          Type = curr_type,
          Seed = k
        )
        
        sim_comparison_tibble <- rbind(sim_comparison_tibble, sim_comparison_tibble_entry)
      }
    }
    
    else {
      
      # Find the current directory
      curr_dir <- paste0(gen_dir, gen_many_seeds_dir, curr_type, "/")
      
      # Read in the tibble containing the data relevant to the current format
      curr_tib <- readRDS(paste0(curr_dir, tib_name))
      
      
      # The relevant index within said tibble
      .curr_index <- which(curr_tib$dataset == d)
      
      # Extract the similarity matrix and expression data
      .curr_sim <- curr_tib$similarity_matrix[[.curr_index]] %>% 
        extract(row_order, row_order)
      
      # .curr_expr <- curr_tib$expression_data[[.curr_index]]
      
      
      .curr_corr <- curr_tib$correlation_matrix[[.curr_index]][row_order, row_order]
      
      # if (curr_type == "500") {
      # Create heatmaps
      ph_sim <- pheatmap(.curr_sim,
                         cluster_rows = F,
                         cluster_cols = F,
                         color = col_pal_sim,
                         breaks = sim_breaks,
                         annotation_row = true_clusters,
                         labels_col = "",
                         silent = T
      )$gtable
      
      # Add to a tibble
      sim_comparison_tibble_entry <- tibble(
        Dataset = d,
        Similarity_matrix = list(.curr_sim),
        Correlation_matrix = list(.curr_corr),
        Expression_data = list(.curr_expr),
        Similarity_matrix_ph = list(ph_sim),
        Correlation_matrix_ph = list(ph_corr),
        Expression_data_ph = list(ph_expr),
        Sim_row_order = list(row_order),
        Expr_col_order = list(expr_col_order),
        Type = curr_type,
        Seed = NA
      )
      
      # Save to big tibble
      sim_comparison_tibble <- rbind(sim_comparison_tibble, sim_comparison_tibble_entry)
    }
  }
}



saveRDS(sim_comparison_tibble, file = "/home/MINTS/sdc56/Desktop/compare_similarity_gen.rds")

# sim_comparison_tibble <- readRDS()

save_dir <- "/home/MINTS/sdc56/Desktop/compare_similarity_gen/"
dir.create(save_dir, showWarnings = F)

sub_dirs <- paste0(save_dir, c("Consensus_500/", "Longer_10/"))

# sim_comparison_tibble <- readRDS("/home/MINTS/sdc56/Desktop/compare_similarity_gen.rds")

primary_type <- "500"

for(dir_ in sub_dirs){
  dir.create(dir_, showWarnings = F)
for (d in datasets) {
  curr_dir <- paste0(dir_, "/", d)
  dir.create(curr_dir, showWarnings = F)

  ph_list <- sim_comparison_tibble$Similarity_matrix_ph[which(sim_comparison_tibble$Dataset == d)]
  types_rel <- sim_comparison_tibble$Type[which(sim_comparison_tibble$Dataset == d)]
  seeds_rel <- sim_comparison_tibble$Seed[which(sim_comparison_tibble$Dataset == d)]
  
  new_order <- c(which(types_rel %in% primary_type), which(! types_rel %in% primary_type))
  ph_list <- ph_list[new_order]
  types_rel <- types_rel[new_order]
  seeds_rel <- seeds_rel[new_order]

  n_ph <- length(ph_list)

  curr_corr_ph <- sim_comparison_tibble$Expression_data_ph[which(sim_comparison_tibble$Dataset == d)][[1]]
  
  if(grepl("Consensus_500", dir_)){
  for (i in 2:n_ph) {
    curr_type <- types_rel[i]
    if (curr_type %in% c("Long", "Longer")) {
      curr_seed <- seeds_rel[i]
      type_save_str <- paste0(curr_type, "_", curr_seed)
      type_string <- paste0(curr_type, " seed ", curr_seed)
    } else {
      type_save_str <- paste0(curr_type, "_many")
      type_string <- paste0("Consensus clustering for ", curr_type, " iterations")
    }

    save_name <- paste0(
      curr_dir,
      "/",
      "500_many_",
      type_save_str,
      ".png"
    )

    main_title <- paste0(
      d,
      "Similarity matrices for consensus clustering for 500 iterations and ",
      type_string
    )

    combine_pheatmaps(list(ph_list[[1]], ph_list[[i]], curr_corr_ph),
      save_name = save_name,
      main = main_title
    )
  }
  } else {
    for (i in 1:(n_ph-1)) {
      curr_type <- types_rel[i]
      if (curr_type %in% c("Long", "Longer")) {
        curr_seed <- seeds_rel[i]
        type_save_str <- paste0(curr_type, "_", curr_seed)
        type_string <- paste0(curr_type, " seed ", curr_seed)
      } else {
        type_save_str <- paste0(curr_type, "_many")
        type_string <- paste0("Consensus clustering for ", curr_type, " iterations")
      }
      
      save_name <- paste0(
        curr_dir,
        "/",
        "longer_10_",
        type_save_str,
        ".png"
      )
      
      main_title <- paste0(
        d,
        "Similarity matrices for longer run (seed 10) and ",
        type_string
      )
      
      combine_pheatmaps(list(ph_list[[n_ph]], ph_list[[i]], curr_corr_ph),
                        save_name = save_name,
                        main = main_title
      )
    }
  }
}

}


# === OLD FOR LOOP =============================================================
# type_defining_order <- "Longer"
# seed_defing_order <- 9
# long_run_seeds <- long_run_seeds[c(which(long_run_seeds == seed_defing_order), which(long_run_seeds != seed_defing_order))]
# for (j in 1:n_types) {
#   curr_type <- types[j]
# 
# 
#   if (curr_type == "Long" | curr_type == "Longer") {
#     # for(curr_type in c("Long", "Longer")){
#     for (k in long_run_seeds) {
#       curr_dir <- paste0(gen_dir, curr_type, "_", k, "/")
# 
#       curr_tib <- readRDS(paste0(curr_dir, tib_name))
# 
#       for (d in datasets) {
#         .curr_index <- which(curr_tib$dataset == d)
# 
#         .curr_sim <- curr_tib$similarity_matrix[[.curr_index]]
#         .curr_expr <- curr_tib$expression_data[[.curr_index]]
# 
#         if(curr_type == type_defining_order & ! is.na(seed_defing_order)){
#           if(k == seed_defing_order){
#             # Find the row order based on the similarity matrix
#             sim_hc <- hclust(dist(.curr_sim))
#             row_order <- sim_hc$order
# 
#             # Find the column order for the expression data
#             expr_hc <- hclust(dist(t(.curr_expr)))
#             expr_col_order <- expr_hc$order
#           }
#         } else {
#           if(! is.na(seed_defing_order)){
#             row_order <- sim_comparison_tibble$Sim_row_order[[which(sim_comparison_tibble$Dataset == d
#                                                                     & sim_comparison_tibble$Type == type_defining_order & sim_comparison_tibble$Seed == seed_defing_order)]]
#             expr_col_order <- sim_comparison_tibble$Expr_col_order[[which(sim_comparison_tibble$Dataset == d
#                                                                           & sim_comparison_tibble$Type == type_defining_order & sim_comparison_tibble$Seed == seed_defing_order)]]
#           } else {
#             row_order <- sim_comparison_tibble$Sim_row_order[[which(sim_comparison_tibble$Dataset == d
#                                                                     & sim_comparison_tibble$Type == type_defining_order)]]
#             expr_col_order <- sim_comparison_tibble$Expr_col_order[[which(sim_comparison_tibble$Dataset == d
#                                                                           & sim_comparison_tibble$Type == type_defining_order)]]
#           }
#         }
# 
#         .curr_sim <- .curr_sim[row_order, row_order]
#         .curr_corr <- curr_tib$correlation_matrix[[.curr_index]][row_order, row_order]
#         .curr_expr <- curr_tib$expression_data[[.curr_index]][row_order, expr_col_order]
# 
#         ph_sim <- pheatmap(.curr_sim,
#                            cluster_rows = F,
#                            cluster_cols = F,
#                            color = col_pal_sim,
#                            breaks = sim_breaks,
#                            annotation_row = true_clusters,
#                            labels_col = "",
#                            silent = T
#         )$gtable
# 
#         ph_corr <- pheatmap(.curr_corr,
#                             cluster_rows = F,
#                             cluster_cols = F,
#                             color = col_pal_corr,
#                             breaks = corr_breaks,
#                             annotation_row = true_clusters,
#                             labels_col = "",
#                             silent = T
#         )$gtable
# 
#         ph_expr <- pheatmap(.curr_expr,
#                             cluster_rows = F,
#                             cluster_cols = F,
#                             color = col_pal_expr,
#                             breaks = expr_breaks,
#                             annotation_row = true_clusters,
#                             labels_col = "",
#                             silent = T
#         )$gtable
# 
#         # Add to a tibble
#         sim_comparison_tibble_entry <- tibble(
#           Dataset = d,
#           Similarity_matrix = list(.curr_sim),
#           Correlation_matrix = list(.curr_corr),
#           Expression_data = list(.curr_expr),
#           Similarity_matrix_ph = list(ph_sim),
#           Correlation_matrix_ph = list(ph_corr),
#           Expression_data_ph = list(ph_expr),
#           Sim_row_order = list(row_order),
#           Expr_col_order = list(expr_col_order),
#           Type = curr_type,
#           Seed = k
#         )
# 
#         sim_comparison_tibble <- rbind(sim_comparison_tibble, sim_comparison_tibble_entry)
#       }
#     }
#   }
# 
#   else {
# 
#     # Find the current directory
#     curr_dir <- paste0(gen_dir, gen_many_seeds_dir, curr_type, "/")
# 
#     # Read in the tibble containing the data relevant to the current format
#     curr_tib <- readRDS(paste0(curr_dir, tib_name))
# 
#     for (d in datasets) {
# 
#       # The relevant index within said tibble
#       .curr_index <- which(curr_tib$dataset == d)
# 
#       # Extract the similarity matrix and expression data
#       .curr_sim <- curr_tib$similarity_matrix[[.curr_index]]
#       .curr_expr <- curr_tib$expression_data[[.curr_index]]
# 
# 
#       if (curr_type == type_defining_order) {
#         # Find the row order based on the similarity matrix
#         sim_hc <- hclust(dist(.curr_sim))
#         row_order <- sim_hc$order
# 
#         # Find the column order for the expression data
#         expr_hc <- hclust(dist(t(.curr_expr)))
#         expr_col_order <- expr_hc$order
#       } else {
#         if(! is.na(seed_defing_order)){
#           row_order <- sim_comparison_tibble$Sim_row_order[[which(sim_comparison_tibble$Dataset == d
#                                                                   & sim_comparison_tibble$Type == type_defining_order & sim_comparison_tibble$Seed == seed_defing_order)]]
#           expr_col_order <- sim_comparison_tibble$Expr_col_order[[which(sim_comparison_tibble$Dataset == d
#                                                                         & sim_comparison_tibble$Type == type_defining_order & sim_comparison_tibble$Seed == seed_defing_order)]]
#         } else {
#           row_order <- sim_comparison_tibble$Sim_row_order[[which(sim_comparison_tibble$Dataset == d
#                                                                   & sim_comparison_tibble$Type == type_defining_order )]]
#           expr_col_order <- sim_comparison_tibble$Expr_col_order[[which(sim_comparison_tibble$Dataset == d
#                                                                         & sim_comparison_tibble$Type == type_defining_order)]]
#         }
#       }
# 
#         # Rearrange data
#         .curr_sim <- .curr_sim[row_order, row_order]
#         .curr_corr <- curr_tib$correlation_matrix[[.curr_index]][row_order, row_order]
#         .curr_expr <- .curr_expr[row_order, expr_col_order]
# 
#         # if (curr_type == "500") {
#         # Create heatmaps
#         ph_sim <- pheatmap(.curr_sim,
#                            cluster_rows = F,
#                            cluster_cols = F,
#                            color = col_pal_sim,
#                            breaks = sim_breaks,
#                            annotation_row = true_clusters,
#                            labels_col = "",
#                            silent = T
#         )$gtable
#         # } else {
#         #
#         #   # Create heatmaps
#         #   ph_sim <- pheatmap(.curr_sim,
#         #     cluster_rows = F,
#         #     cluster_cols = F,
#         #     color = col_pal_sim,
#         #     breaks = sim_breaks,
#         #     labels_col = "",
#         #     silent = T
#         #   )$gtable
#         # }
# 
#         ph_corr <- pheatmap(.curr_corr,
#                             cluster_rows = F,
#                             cluster_cols = F,
#                             color = col_pal_corr,
#                             breaks = corr_breaks,
#                             labels_col = "",
#                             silent = T
#         )$gtable
# 
#         ph_expr <- pheatmap(.curr_expr,
#                             cluster_rows = F,
#                             cluster_cols = F,
#                             color = col_pal_expr,
#                             breaks = expr_breaks,
#                             annotation_row = true_clusters,
#                             labels_col = "",
#                             silent = T
#         )$gtable
# 
#         # Add to a tibble
#         sim_comparison_tibble_entry <- tibble(
#           Dataset = d,
#           Similarity_matrix = list(.curr_sim),
#           Correlation_matrix = list(.curr_corr),
#           Expression_data = list(.curr_expr),
#           Similarity_matrix_ph = list(ph_sim),
#           Correlation_matrix_ph = list(ph_corr),
#           Expression_data_ph = list(ph_expr),
#           Sim_row_order = list(row_order),
#           Expr_col_order = list(expr_col_order),
#           Type = curr_type,
#           Seed = NA
#         )
# 
#         # Save to big tibble
#         sim_comparison_tibble <- rbind(sim_comparison_tibble, sim_comparison_tibble_entry)
#       }
#     }
#   }
#   
# 
# save_dir <- "/home/MINTS/sdc56/Desktop/compare_similarity_gen/"
# dir.create(save_dir, showWarnings = F)
# 
# 
# 
# # sim_comparison_tibble <- readRDS("/home/MINTS/sdc56/Desktop/compare_similarity_gen.rds")
# 
# primary_type <- "Longer"
# primary_seed <- 9
# primary_type_save_str <- "Longer_9"
# # type_save_str <- paste0(curr_type, "_", curr_seed)
# primary_type_string <- paste0(primary_type, " seed ", primary_seed)
# 
# sub_dirs <- paste0(save_dir, primary_type)
# 
# for(dir_ in sub_dirs){
#   dir.create(dir_, showWarnings = F)
#   for (d in datasets) {
#     curr_dir <- paste0(dir_, "/", d)
#     dir.create(curr_dir, showWarnings = F)
#     
#     ph_list <- sim_comparison_tibble$Similarity_matrix_ph[which(sim_comparison_tibble$Dataset == d)]
#     types_rel <- sim_comparison_tibble$Type[which(sim_comparison_tibble$Dataset == d)]
#     seeds_rel <- sim_comparison_tibble$Seed[which(sim_comparison_tibble$Dataset == d)]
#     
#     new_order <- c(which(types_rel %in% primary_type), which(! types_rel %in% primary_type))
#     ph_list <- ph_list[new_order]
#     types_rel <- types_rel[new_order]
#     seeds_rel <- seeds_rel[new_order]
#     
#     n_ph <- length(ph_list)
#     
#     curr_corr_ph <- sim_comparison_tibble$Expression_data_ph[which(sim_comparison_tibble$Dataset == d)][[1]]
#     
#     
#     
#     # if(grepl("Consensus_500", dir_)){
#       for (i in 2:n_ph) {
#         curr_type <- types_rel[i]
#         if (curr_type %in% c("Long", "Longer")) {
#           curr_seed <- seeds_rel[i]
#           type_save_str <- paste0(curr_type, "_", curr_seed)
#           type_string <- paste0(curr_type, " seed ", curr_seed)
#         } else {
#           type_save_str <- paste0(curr_type, "_many")
#           type_string <- paste0("Consensus clustering for ", curr_type, " iterations")
#         }
#         
#         save_name <- paste0(
#           curr_dir,
#           "/",
#           primary_type_save_str, 
#           "_",
#           type_save_str,
#           ".png"
#         )
#         
#         main_title <- paste0(
#           d, "Similarity matrices for ", primary_type_string, " and ", type_string
#         )
#         
#         combine_pheatmaps(list(ph_list[[1]], ph_list[[i]], curr_corr_ph),
#                           save_name = save_name,
#                           main = main_title
#         )
#       }
#   }
#   
# }
