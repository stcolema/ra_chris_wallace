

library(mcclust)
library(magrittr)
library(ggplot2)
library(ggforce)


# For tibbles
library(tibble) # for dataframe of lists

# For data wrangling
library(dplyr)

# Heatmapping
library(pheatmap) # install.packages("pheatmap", dep = T)

# Colour palettes
library(RColorBrewer)

# To load true clustering
library(rlist)


function_dir <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/Analysis_script_functions/"

function_scripts <- c(
  "plot_rand_index.R",
  "plot_comparison_expression_clustering.R"
)

for (f in paste0(function_dir, function_scripts)) {
  source(f)
}

# For controlling ordering of boxplots
add_level <- function(x, new_level){
  if(is.factor(x)) return(factor(x, levels=c(levels(x), new_level)))
  return(x)
}

# === Main script ==============================================================

# Set ggplot2 theme
theme_set(theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust=1)))

true_clusters <- readRDS("~/Desktop/ra_chris_wallace/Yeast/true_clusters_for_yeast_data.rds")

my_tib <- readRDS("~/Desktop/Yeast_3_datasets/MDI_output/Consensus_500/compare_tibble.rds")

ann_row <- data.frame(Cluster = as.factor(true_clusters[[1]]))


col_pal <- colorRampPalette(c("#FF9900", "white", "#146EB4"))(100)
my_breaks <- define_breaks(col_pal)

cor_pal <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(100)
cor_breaks <- define_breaks(cor_pal)

curr_cor <- my_tib$correlation_matrix[[1]]
row.names(curr_cor) <- row.names(ann_row)
ann_colours <- list(Cluster = c(1 = ""))

ann_colours <- c("#fbb4ae",
  "#b3cde3",
  "#ccebc5",
  "#decbe4",
  "#fed9a6",
  "#ffffcc",
  "#e5d8bd",
  "#fddaec",
  "#1f78b4")
  # "#f2f2f2")

names(ann_colours) =1:9

ann_colours_dt <- list(Cluster = ann_colours)

# row_order <- match(row.names(ann_row), row.names(my_tib$correlation_matrix[[1]]))
pheatmap(curr_cor,
         cluster_rows = F,
         cluster_cols = F,
         annotation_row = ann_row,
         annotation_colors = ann_colours_dt, 
         color = cor_pal,
         breaks = cor_breaks,
         show_rownames = F,
         show_colnames = F,
         main = "Simulation 1: correlation matrix for dataset 1",
         filename = "~/Desktop/ra_chris_wallace/Notes/Thesis/Images/Gen_data/Case_1/cor_matrix_dataset_1.png")
