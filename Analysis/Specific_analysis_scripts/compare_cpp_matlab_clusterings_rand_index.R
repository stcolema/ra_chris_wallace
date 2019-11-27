

library(tibble)
library(ggplot2)
library(mclust)
library(magrittr)
library(pheatmap)

input_arguments <- function() {
  option_list <- list(

    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--datasets"),
      type = "character",
      default = "CD14.csv CD15.csv CD19.csv CD4.csv CD8.csv IL.csv PLA.csv RE.csv TR.csv",
      help = "command to transpose all [EXT] files in currect directory [default= %default]",
      metavar = "character"
    ),

    # Directory to read from
    optparse::make_option(c("-d", "--dir"),
      type = "character",
      default = ".",
      help = "directory to read files from [default= %default]",
      metavar = "character"
    ),

    # Directory to read from
    optparse::make_option(c("-c", "--cpp_tibble"),
      type = "character",
      default = NA,
      help = "file path to tibble containing data about the C++ MDI [default= %default]",
      metavar = "character"
    ),

    # Directory to read from
    optparse::make_option(c("-m", "--matlab_tibble"),
      type = "character",
      default = NA,
      help = "file path to tibble containing data about the MATLAB MDI [default= %default]",
      metavar = "character"
    ),

    # Directory to read from
    optparse::make_option(c("-t", "--title"),
      type = "character",
      default = "Distribution of adjusted rand index between C++ and MATLAB MDI",
      help = "Main title for output plot [default= %default]",
      metavar = "character"
    ),

    # Directory to read from
    optparse::make_option(c("-s", "--subtitle"),
      type = "character",
      default = NA,
      help = "Subtitle for output plot [default= %default]",
      metavar = "character"
    ),

    # Directory to read from
    optparse::make_option(c("--save_name"),
      type = "character",
      default = "./comparison_rand_index.png",
      help = "Save name for output plot [default= %default]",
      metavar = "character"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}


# === Main script ==============================================================

# Read input arguments
args <- input_arguments()

# Extract file names
many_seeds_tbl_filename <- args$cpp_tibble
matlab_tbl_filename <- args$matlab_tibble

# Extract plot titles
plot_title <- args$title
plot_subtitle <- args$subtitle

# Read in tibbles
many_seeds_tbl <- readRDS(many_seeds_tbl_filename)
matlab_tbl <- readRDS(matlab_tbl_filename)

# Output file name
save_name <- args$save_name

# The datasets
datasets <- many_seeds_tbl$dataset

# Various numbers for sizing and looping
n_datasets <- length(datasets)
n_seeds <- nrow(many_seeds_tbl$mdi_allocation[[1]])
n_comparisons <- n_datasets * n_seeds

# Create the vector of dataset names
dataset_vec <- lapply(datasets, rep, n_seeds) %>% unlist()

# Create the vector of seeds
seed_vec <- rep(1:n_seeds, n_datasets)

# Data frame to hold adjusted rand index for a given seed of C++ MDI and a given
# dataset compared to the MATLAB clustering for the same dataset
rand_df <- data.frame(
  Dataset = dataset_vec,
  Seed = seed_vec,
  Rand_index = rep(0, n_comparisons)
)
# Iterate and make rand index
for (i in 1:n_comparisons) {

  # The current dataset
  curr_dataset <- rand_df$Dataset[i]

  # Current seed for C++
  curr_seed <- rand_df$Seed[i]

  # The entry in the MATLAB tibble
  matlab_tibble_index <- matlab_tbl$dataset == curr_dataset

  # The MATLAB clustering
  curr_matlab_clustering <- matlab_tbl$pred_allocation %>%
    magrittr::extract(matlab_tibble_index) %>%
    magrittr::extract2(1)
    
  # curr_matlab_clustering <- 
  #        matlab_tbl$similarity_matrix %>% 
  #        magrittr::extract(matlab_tibble_index) %>%
  #        magrittr::extract2(1) %>% 
  #        round() %>% 
  #        mcclust::Simtocl() 

  # The Many seeds tibble row entry
  mason_tibble_index_outer <- many_seeds_tbl$dataset == curr_dataset

  # Within the MDI allocations the row to compare to the MATLAB
  mason_tibble_index_inner <- curr_seed

  # Extract the MDI custering for the approrpiate dataset and seed
  curr_mason_clustering <- many_seeds_tbl$mdi_allocation %>%
    magrittr::extract(mason_tibble_index_outer) %>%
    magrittr::extract2(1) %>%
    magrittr::extract(mason_tibble_index_inner, ) %>%
    unlist()
  
  # Compute the addjusted rand index
  curr_rand_index <- adjustedRandIndex(curr_matlab_clustering, curr_mason_clustering)

  # Record
  rand_df$Rand_index[i] <- curr_rand_index
}

# Plot
ggplot(data = rand_df, aes(x = Rand_index)) +
  geom_histogram() +
  facet_wrap(~Dataset) +
  labs(
    title = plot_title,
    subtitle = plot_subtitle,
    x = "Adjusted rand index",
    y = "Count"
  )

ggsave(save_name)
