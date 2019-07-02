





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

# For str_replace
library(stringr)

# === Main script ==============================================================

# Directroy containing DATA, MDI output and other stuff
gen_dir <- "/home/MINTS/sdc56/Desktop/Gen_data_output/"
use_predefined_row_order <- F
type_defining_order <- "Long_8"
n_genes <- 400


mdi_main_dir <- paste0(gen_dir, "MDI_output")
mdi_dirs <- list.dirs(mdi_main_dir, recursive = F)
types <- list.dirs(mdi_main_dir, recursive = F, full.names = F)

tib_name <- paste0(mdi_dirs, "/compare_tibble.rds")
all_tibbles <- lapply(tib_name, readRDS) %>% 
  set_names(types)

# types <- c("Longer", "500", "1000", "5000", "10000", "Long")
# long_run_seeds <- c(9, 1:8, 10)

# n_weights <- length(weights)
n_types <- length(types)

datasets <- c("MDItestdata1", "MDItestdata2", "MDItestdata3")

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
  Type = character()
  )







# === Heatmap inputs ===========================================================

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

# Annotation labels based on generating method
true_clusters <- data.frame(
  Labels = as.factor(c(rep(1, 25), rep(2, 50), rep(3, 75), rep(4, 100), rep(5, 150)))
) %>% set_rownames(paste0("Gene_", 1:n_genes))

if(use_predefined_row_order){
  row_order <- order(true_clusters$Labels)
}

# === Record similarity matrices ===============================================

# first_type <- which(types %in% type_defining_order)
# types <- types[c(which(types %in% type_defining_order), which(! types %in% type_defining_order))]
# seed_defing_order <- NA 

main_tib <- all_tibbles[[type_defining_order]]
for (d in datasets) {
  # The relevant index within said tibble
  
  
  .curr_index <- which(datasets == d)
  

  if(! use_predefined_row_order){
    sim_hc <- main_tib$similarity_matrix[[.curr_index]] %>% 
      t() %>% 
      dist() %>% 
      hclust()
    row_order <- sim_hc$order
  }
  
  # .curr_expr <- curr_tib$expression_data[[.curr_index]]
  .curr_expr <- main_tib$expression_data[[.curr_index]]
  
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
    curr_tib <- all_tibbles[[curr_type]]

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
      Type = curr_type
    )
    
    sim_comparison_tibble <- rbind(sim_comparison_tibble, sim_comparison_tibble_entry)
  }
}

# === Save comparison plots ====================================================


save_dir <- paste0(gen_dir, "Compare_similarity_matrices")
dir.create(save_dir, showWarnings = F)

sub_dir <- paste0(save_dir, "/", type_defining_order)

re_ordered_types <- types[
  c(which(types == type_defining_order), which(types != type_defining_order))
]

if(grepl("Long", type_defining_order)){
  type_defining_order_title_str <- str_replace(type_defining_order, "_", " run: seed ")
} else {
  type_defining_order_title_str <- str_replace(type_defining_order, "Many_seeds_", "Consensus ")
}

# for(sub_dir in sub_dirs){
dir.create(sub_dir, showWarnings = F)
for (d in datasets) {
  curr_dir <- paste0(sub_dir, "/", d)
  dir.create(curr_dir, showWarnings = F)
  
  ph_list <- sim_comparison_tibble$Similarity_matrix_ph[which(sim_comparison_tibble$Dataset == d)]
  types_rel <- sim_comparison_tibble$Type[which(sim_comparison_tibble$Dataset == d)]

  new_order <- c(which(types_rel == type_defining_order), which(types_rel != type_defining_order))
  ph_list <- ph_list[new_order]
  types_rel <- types_rel[new_order]

  n_ph <- length(ph_list)
  
  curr_corr_ph <- sim_comparison_tibble$Expression_data_ph[which(sim_comparison_tibble$Dataset == d)][[1]]
  
# if(grepl("Consensus_500", sub_dir)){
  for (i in 2:n_ph) {
    curr_type <- types_rel[i]
    if (grepl("Long", curr_type)) {
      type_string <- str_replace(curr_type, "_", " run: seed ")
      type_save_str <- paste0(curr_type)
    } else {
      type_save_str <- paste0(curr_type, "_many")
      type_string <- str_replace(curr_type, "Many_seeds_", "Consensus ")
    }
    
    save_name <- paste0(
      curr_dir,
      "/",
      type_defining_order,
      "_",
      type_save_str,
      ".png"
    )
    
    main_title <- paste0(
      d,
      ": Similarity matrices for ",
      type_defining_order_title_str,
      " and ",
      type_string
    )
    
    combine_pheatmaps(list(ph_list[[1]], ph_list[[i]], curr_corr_ph),
                      save_name = save_name,
                      main = main_title
    )
  }
}
# }
