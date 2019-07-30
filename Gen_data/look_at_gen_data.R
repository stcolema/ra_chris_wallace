
# Script inspecting the generated data and it's labelling

library(data.table)
library(magrittr)
library(pheatmap)
library(RColorBrewer)

function_dir <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/Analysis_script_functions/"

function_scripts <- c(
  "plot_comparison_expression_clustering.R"
)

for (f in paste0(function_dir, function_scripts)) {
  source(f)
}

# Read in the data
d1 <- read.csv("~/Desktop/Gen_data/mk_0/new_data.csv", header = T, row.names = "V1")
d2 <- read.csv("~/Desktop/Gen_data/mk_1/new_data.csv", header = T, row.names = "V1")
d3 <- read.csv("~/Desktop/Gen_data/mk_2/new_data.csv", header = T, row.names = "V1")

# The true labels based on how the data is generated
true_labels <- c(
  rep(1, 25),
  rep(2, 50),
  rep(3, 75),
  rep(4, 100),
  rep(5, 150)
)


# Annotation data
annotation_row <- data.frame(Cluster = as.factor(true_labels)) %>% # , levels = 1:5)) %>%
  set_rownames(row.names(d1))

# Annotation colours
var_1 <- brewer.pal(length(unique(true_labels)), "BrBG")

names(var_1) <- paste0("Cluster", 1:5)
names(var_1) <- 1:5
col_pal <- list(Class = var_1)

ann_colours <- c("#fbb4ae",
  "#b3cde3",
  "#ccebc5",
  "#decbe4",
  "#fed9a6") %>% 
  set_names(1:5)

ann_colours_df <- list(Cluster = ann_colours)

gen_file_name <- "~/Desktop/ra_chris_wallace/Notes/Thesis/Images/Gen_data/Case_2/dataset_"

# === Heatmaps =================================================================

# These are not saved. If called from terminal these heatmaps will be saved to 
# Rplots.pdf together
pheatmap(d1,
  annotation_row = annotation_row,
  cluster_rows = F,
  cluster_cols = F,
  filename = paste0(gen_file_name, 1, ".png"),
  annotation_colors = ann_colours_df,
  main = "Simulation 2: Dataset 1",
  show_rownames = F,
  show_colnames = F
)

pheatmap(d2,
  annotation_row = annotation_row,
  filename = paste0(gen_file_name, 2, ".png"),
  cluster_rows = F,
  cluster_cols = F,
  annotation_colors = ann_colours_df,
  main = "Simulation 2: Dataset 2",
  show_rownames = F,
  show_colnames = F
)
pheatmap(d3,
  annotation_row = annotation_row, 
  filename = paste0(gen_file_name, 3, ".png"),
  cluster_rows = F,
  cluster_cols = F,
  annotation_colors = ann_colours_df,
  main = "Simulation 2: Dataset 3",
  show_rownames = F,
  show_colnames = F
)

# These are not saved. If called from terminal these heatmaps will be saved to 
# Rplots.pdf together
pheatmap(d1,
         annotation_row = annotation_row,
         cluster_cols = F,
         filename = paste0(gen_file_name, 1, "_clustered_rows.png"),
         cluster_rows = T,
         main = "Simulation 2: Dataset 1 (clustered rows)",
         show_rownames = F,
         show_colnames = F
)

pheatmap(d2,
         annotation_row = annotation_row,
         filename = paste0(gen_file_name, 2, "_clustered_rows.png"),
         cluster_rows = T,
         cluster_cols = F,
         annotation_colors = ann_colours_df,
         main = "Simulation 2: Dataset 2 (clustered rows)",
         show_rownames = F,
         show_colnames = F
)
pheatmap(d3,
         annotation_row = annotation_row, 
         filename = paste0(gen_file_name, 3, "_clustered_rows.png"),
         cluster_rows = T,
         cluster_cols = F,
         annotation_colors = ann_colours_df,
         main = "Simulation 2: Dataset 3 (clustered rows)",
         show_rownames = F,
         show_colnames = F
)

ph3 <- pheatmap(d3,
                annotation_row = annotation_row, 
                # filename = paste0(gen_file_name, 3, ".png"),
                cluster_rows = F,
                cluster_cols = F,
                annotation_colors = ann_colours_df,
                # main = "Simulation 2: Dataset 3",
                show_rownames = F,
                show_colnames = F,
                silent = T
)$gtable

row_order <- hclust(dist(d3))$order

ph3_ordered <- pheatmap(d3[row_order,],
         annotation_row = annotation_row, 
         # filename = paste0(gen_file_name, 3, "_clustered_rows.png"),
         cluster_rows = F,
         cluster_cols = F,
         annotation_colors = ann_colours_df,
         # main = "Simulation 2: Dataset 3 (clustered rows)",
         show_rownames = F,
         show_colnames = F,
         silent = T
)$gtable

combine_pheatmaps(list(ph3, ph3_ordered),
                  save_name = paste0(gen_file_name, 3, "_comp_clustered_unclustered.png"),
                  main = "Simulation 2: Dataset 3")




pheatmap(d1,
  annotation_row = annotation_row, cluster_cols = F, cluster_rows = F,
  # filename = paste0(gen_file_name, 1, ".pdf"),
  annotation_colors = col_pal,
  main = "Dataset 1"
)
pheatmap(d2,
  annotation_row = annotation_row, cluster_cols = F, cluster_rows = F,
  # filename = paste0(gen_file_name, 2, ".pdf"),
  annotation_colors = col_pal,
  main = "Dataset 2"
)
pheatmap(d3,
  annotation_row = annotation_row, cluster_cols = F, cluster_rows = F,
  # filename = paste0(gen_file_name, 3, ".pdf"),
  annotation_colors = col_pal,
  main = "Dataset 3"
)

# === Saving heatmaps ==========================================================

gen_file_name <- "~/Desktop/pheatmap_data_"

# Clustered rows
pheatmap(d1,
  annotation_row = annotation_row,
  cluster_cols = F,
  filename = paste0(gen_file_name, 1, "_clustered_rows.pdf"),
  annotation_colors = col_pal,
  main = "Dataset 1"
)
pheatmap(d2,
  annotation_row = annotation_row, cluster_cols = F,
  filename = paste0(gen_file_name, 2, "_clustered_rows.pdf"),
  annotation_colors = col_pal,
  main = "Dataset 2"
)
pheatmap(d3,
  annotation_row = annotation_row, cluster_cols = F,
  filename = paste0(gen_file_name, 3, "_clustered_rows.pdf"),
  annotation_colors = col_pal,
  main = "Dataset 3"
)

# Unclustered rows
pheatmap(d1,
  annotation_row = annotation_row, cluster_cols = F, cluster_rows = F,
  filename = paste0(gen_file_name, 1, ".pdf"),
  annotation_colors = col_pal,
  main = "Dataset 1"
)
pheatmap(d2,
  annotation_row = annotation_row, cluster_cols = F, cluster_rows = F,
  filename = paste0(gen_file_name, 2, ".pdf"),
  annotation_colors = col_pal,
  main = "Dataset 2"
)
pheatmap(d3,
  annotation_row = annotation_row, cluster_cols = F, cluster_rows = F,
  filename = paste0(gen_file_name, 3, ".pdf"),
  annotation_colors = col_pal,
  main = "Dataset 3"
)
