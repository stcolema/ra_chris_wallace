
# Script inspecting the generated data and it's labelling

library(data.table)
library(magrittr)
library(pheatmap)
library(RColorBrewer)

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
annotation_row <- data.frame(Class = as.factor(true_labels)) %>% # , levels = 1:5)) %>%
  set_rownames(row.names(d1))

# Annotation colours
var_1 <- brewer.pal(length(unique(true_labels)), "BrBG")

names(var_1) <- paste0("C", 1:5)
names(var_1) <- 1:5
col_pal <- list(Class = var_1)

# === Heatmaps =================================================================

# These are not saved. If called from terminal these heatmaps will be saved to 
# Rplots.pdf together
pheatmap(d1,
  annotation_row = annotation_row,
  cluster_cols = F,
  # filename = paste0(gen_file_name, 1, "_clustered_rows.pdf"),
  annotation_colors = col_pal,
  main = "Dataset 1"
)
pheatmap(d2,
  annotation_row = annotation_row, cluster_cols = F,
  # filename = paste0(gen_file_name, 2, "_clustered_rows.pdf"),
  annotation_colors = col_pal,
  main = "Dataset 2"
)
pheatmap(d3,
  annotation_row = annotation_row, cluster_cols = F,
  # filename = paste0(gen_file_name, 3, "_clustered_rows.pdf"),
  annotation_colors = col_pal,
  main = "Dataset 3"
)

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
