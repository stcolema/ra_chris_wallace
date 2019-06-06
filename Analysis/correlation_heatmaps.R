library(pheatmap)
library(RColorBrewer)
library(data.table)
library(magrittr)
library(optparse)

input_arguments <- function() {
  option_list <- list(

    # Directory to read from
    optparse::make_option(c("-d", "--dir"),
      type = "character",
      default = ".",
      help = "directory to read files from [default= %default]",
      metavar = "character"
    ),

    optparse::make_option(c("-w", "--write_dir"),
      type = "character",
      default = NULL,
      help = "directory to write files to [defaults to --dir]",
      metavar = "character"
    ),

    optparse::make_option(c("-p", "--plot_type"),
      type = "character",
      default = ".png",
      help = "file type of heatmaps (one of '.pdf' or '.png') [default= %default]",
      metavar = "character"
    ),

    optparse::make_option(c("-k", "--probe_key"),
      type = "character",
      default = "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/probe_key_unique_names.csv",
      help = "file to read gene names from [default= %default]",
      metavar = "character"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

paul_pheatmap <- function(datasetToPlot, save_name = NA) {
  paletteLength <- 50
  myColor <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(paletteLength)
  myBreaks <- c(
    seq(min(datasetToPlot), 0, length.out = ceiling(paletteLength / 2) + 1),
    seq(max(datasetToPlot) / paletteLength, max(datasetToPlot), length.out = floor(paletteLength / 2))
  )

  ph <- pheatmap(datasetToPlot,
    cluster_rows = T,
    cluster_cols = T,
    color = myColor,
    breaks = myBreaks,
    filename = save_name
  )
  ph
}

save_heatmaps <- function(data_files, save_names, probe_key) {
  n_files <- length(data_files)

  # Colours and breaks for correlation heatmap
  palette_length <- 50
  my_colour <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(palette_length)
  my_breaks <- c(
    seq(-1, 0, length.out = ceiling(palette_length / 2) + 1),
    seq(1 / palette_length, 1, length.out = floor(palette_length / 2))
  )

  for (i in 1:n_files) {

    # Find the current file and relevant save names
    curr_expr_file <- data_files[[i]]
    curr_save_names <- save_names[[i]]

    # Read in the data
    raw_data <- fread(curr_expr_file, header = T)

    # Remove the V1 column containing probe IDs and set the row names
    print(ncol(raw_data))
    
    my_df <- raw_data[, -c(1)] %>%
      set_rownames(unlist(raw_data[, 1]))

    # Find the corresponding gene_names
    gene_names <- probe_key$Unique_gene_name[match(row.names(my_df), probe_key$ProbeID)]

    print(length(gene_names))
    
    # Set the row names as gene names
    row.names(my_df) <- gene_names

    # Create the correlation matrix
    cor_mat <- cor(t(my_df))

    # Heatmap
    ph_cor <- pheatmap(cor_mat,
      color = my_colour,
      breaks = my_breaks,
      filename = curr_save_names[[1]]
    )


    paul_pheatmap(my_df, save_name = curr_save_names[[2]])
  }
}

# Read in arguments
args <- input_arguments()
data_dir <- args$dir
write_dir <- args$write_dir
plot_type <- args$plot_type
probe_key <- fread(args$probe_key)

acceptable_plot_types <- c(".pdf", ".png")

if (!plot_type %in% acceptable_plot_types) {
  stop("--plot_type must be one of '.png' or '.pdf'.")
}

# If NULL values for write directory, write to same location as read directory
if (is.null(write_dir)) {
  write_dir <- data_dir
}

# If it does not already exist, create the write directory
dir.create(write_dir, showWarnings = FALSE)

# List the files of interest
curr_names <- list.files(path = data_dir, full.names = T, include.dirs = F) %>%
  grep("csv", ., value = TRUE)

# Reduced file name
file_names <- basename(tools::file_path_sans_ext(curr_names))


# Create save locations for the different heatmaps
gen_save_dir <- paste(write_dir, file_names)

for (d in gen_save_dir) {
  dir.create(d, showWarnings = FALSE)
}

heatmap_names <- c(
  "Correlation",
  "Expression"
)

full_save_names <- gen_save_dir %>%
  lapply(paste0, "/", heatmap_names, plot_type)

save_heatmaps(curr_names, full_save_names, probe_key)
