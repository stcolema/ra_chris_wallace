

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
  pheatmap(matlab_sim[row_order, row_order],
    cluster_rows = F,
    cluster_cols = F,
    main = paste0(dataset, ": PSM for MATLAB"),
    filename = ph_filenames[2],
    cellheight = 3.0
  )

  # Heatmap of expression data
  pheatmap(expression_data[row_order, ],
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

