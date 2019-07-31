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

# imnputs for heatmaps
col_pal <- colorRampPalette(c("#FF9900", "white", "#146EB4"))(100)
my_breaks <- define_breaks(col_pal)

cor_pal <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(100)
cor_breaks <- define_breaks(cor_pal)


tib_consensus <- readRDS("~/Desktop/Set_250/Specific_datasets/CD14/compare_tibble.rds")
tib_bayes_1 <- readRDS("~/Desktop/Mixture_models/CD14/out_seed_1/compare_tibble.rds") 
tib_bayes_2 <- readRDS("~/Desktop/Mixture_models/CD14/out_seed_2/compare_tibble.rds") 
tib_bayes_3 <- readRDS("~/Desktop/Mixture_models/CD14/out_seed_3/compare_tibble.rds") 


con_sim <- tib_consensus$similarity_matrix[[1]]
bayes_sim_1 <- tib_bayes_1$similarity_matrix[[1]]
bayes_sim_2 <- tib_bayes_2$similarity_matrix[[1]]
bayes_sim_3 <- tib_bayes_3$similarity_matrix[[1]]

row_order <- hclust(dist(con_sim))$order
ph_1 <- pheatmap(con_sim[row_order, row_order],
                 color = col_pal, 
                 breaks = my_breaks,
                 cluster_rows = F,
                 cluster_cols = F)

ph_2 <- pheatmap(bayes_sim_1[row_order, row_order],
                 color = col_pal, 
                 breaks = my_breaks,
                 cluster_rows = F,
                 cluster_cols = F)


# Bayes ordering 
b_row_order <- hclust(dist(bayes_sim_1))$order



ph_bayes_1 <- pheatmap(bayes_sim_1[b_row_order, b_row_order],
                 color = col_pal, 
                 breaks = my_breaks,
                 cluster_rows = F,
                 cluster_cols = F,
                 show_rownames = F,
                 show_colnames = F,
                 silent = T)$gtable

ph_bayes_2 <- pheatmap(bayes_sim_2[b_row_order, b_row_order],
                 color = col_pal, 
                 breaks = my_breaks,
                 cluster_rows = F,
                 cluster_cols = F,
                 show_rownames = F,
                 show_colnames = F,
                 silent = T)$gtable

ph_bayes_3 <- pheatmap(bayes_sim_3[b_row_order, b_row_order],
                 color = col_pal, 
                 breaks = my_breaks,
                 cluster_rows = F,
                 cluster_cols = F,
                 show_rownames = F,
                 show_colnames = F,
                 silent = T)$gtable



combine_pheatmaps(list(ph_bayes_1, ph_bayes_2, ph_bayes_3),
                              save_name = "~/Desktop/ra_chris_wallace/Notes/Thesis/Images/Biology_data/Set_250/Bayesian_mixture_models/CD14_comparison_across_seeds.png",
                              main = "CEDAR 1: CD14: Comparison of Bayesian mixture models",
                              font_size = 20)
