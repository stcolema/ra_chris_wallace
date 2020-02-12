#!/usr/bin/env Rscript

# For tibbles
library(tibble) # for dataframe of lists

# For data wrangling
library(dplyr)

# Heatmapping
library(pheatmap) # install.packages("pheatmap", dep = T)

# Colour palettes
library(RColorBrewer)

# For command line arguments
library(optparse) # install.packages("optparse")

# Load data.table to access fread and fwrite functions
library(data.table) # install.packages("data.table", dep = T)

# Load magrittr for the pipe %>%
library(magrittr)


input_arguments <- function() {
  option_list <- list(

    # # Number of genes to generate
    # optparse::make_option(c("-n", "--n_rows"),
    #                       type = "integer",
    #                       default = 400L,
    #                       help = "Number of genes (rows) to generate [default= %default]",
    #                       metavar = "integer"
    # ),

    # Number of people to generate
    optparse::make_option(c("-p", "--n_cols"),
      type = "integer",
      default = 100L,
      help = "Number of people (columns) to generate [default= %default]",
      metavar = "integer"
    ),

    # # Number of genes to generate
    # optparse::make_option(c("-n"),
    #                       type = "integer",
    #                       default = 400L,
    #                       help = "Number of genes to generate [default= %default]",
    #                       metavar = "integer"
    # ),

    # Number of genes to generate
    optparse::make_option(c("-m", "--mu"),
      type = "numeric",
      default = 0,
      help = "Mean of data generated to peturb data [default= %default]",
      metavar = "integer"
    ),

    # Number of genes to generate
    optparse::make_option(c("-s", "--sd"),
      type = "numeric",
      default = 0.1,
      help = "Number of genes to generate [default= %default]",
      metavar = "integer"
    ),

    # Number of genes to generate in each cluster
    optparse::make_option(c("--n_clust"),
      type = "character",
      default = "25 50 75 100 150",
      help = "Number of genes to generate in each cluster [default= %default]",
      metavar = "character"
    ),

    # Means of generated clusters
    optparse::make_option(c("--means"),
      type = "character",
      default = "1 2 4 6 8",
      help = "Means of generated clusters [default= %default]",
      metavar = "character"
    ),

    # Instuction to save pheatmap of data generated
    optparse::make_option(c("--plot_heatmap"),
      type = "logical",
      default = TRUE,
      help = "Instuction to save pheatmap of data generated [default= %default]",
      metavar = "logical"
    ),

    # Filename of pheatmap (if saved)
    optparse::make_option(c("--ph_save_name"),
      type = "character",
      default = "data_pheatmap.png",
      help = "Filename of pheatmap (if saved) [default= %default]",
      metavar = "character"
    ),

    # Instuction to save pheatmap of data generated
    optparse::make_option(c("--plot_cor_heatmap"),
      type = "logical",
      default = TRUE,
      help = "Instuction to save pheatmap of correlation within generated data [default= %default]",
      metavar = "logical"
    ),

    # Filename of pheatmap (if saved)
    optparse::make_option(c("--cor_ph_save_name"),
      type = "character",
      default = "cor_pheatmap.png",
      help = "Filename of correlation pheatmap (if saved) [default= %default]",
      metavar = "character"
    ),

    # Directory to save data and heatmap (if saved)
    optparse::make_option(c("-d", "--dir"),
      type = "character",
      default = "./",
      help = "Directory to save data and heatmap (if saved) [default= %default]",
      metavar = "character"
    ),

    # Filename for data
    optparse::make_option(c("-f", "--file"),
      type = "character",
      default = "new_data.csv",
      help = "Filename for data [default= %default]",
      metavar = "character"
    ),

    # Random seed to use
    optparse::make_option(c("--seed"),
      type = "numeric",
      default = 1,
      help = "Random seed [default= %default]",
      metavar = "integer"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}


generate_peturbed_data <- function(means,
                                   n_clust,
                                   mu_peturbation,
                                   sd_peturbation,
                                   p) {

  # Basic data of 5 clusters
  univariateData <- c(
    rnorm(n_clust[1], means[1]),
    rnorm(n_clust[2], means[2]),
    rnorm(n_clust[3], means[3]),
    rnorm(n_clust[4], means[4]),
    rnorm(n_clust[5], means[5])
  )

  # Declare generated data
  gen_data <- c()
  for (i in 1:p) {
    gen_data <- rbind(gen_data, univariateData + rnorm(n, mu_peturbation, sd_peturbation))
  }

  # Transpose
  t(gen_data)
}

# === Receive command line arguments ===========================================

args <- input_arguments()

# The random seed to enable reproducibility
seed <- args$seed

# The number of columns to generate
p <- args$n_cols

# The parameters to define the peturbation to ensure less defined clusters
# The mean peturbation
mu_peturbation <- args$mu

# The stnadard deviation of the peturbation
sd_peturbation <- args$sd

# The number of clusters to generate in each cluster
n_clust <- args$n_clust %>%
  strsplit(" ") %>%
  unlist() %>%
  as.numeric()

# The means of the clusters
means <- args$means %>%
  strsplit(" ") %>%
  unlist() %>%
  as.numeric()

# Directory to save to
loc_dir <- args$dir

# Filename of generated data
filename <- args$file

# Save name for heatmap
ph_save_name <- args$ph_save_name

# Instruction to make a heatmap
do_heatplot <- args$plot_heatmap

# Save name for correlation heatmap
cor_ph_save_name <- args$cor_ph_save_name

# Instruction to make a correlation heatmap
do_cor_heatplot <- args$plot_cor_heatmap

# === Main script ==============================================================

set.seed(seed)
n <- sum(n_clust)

# Generate data
new_data <- generate_peturbed_data(
  means,
  n_clust,
  mu_peturbation,
  sd_peturbation,
  p
)

# Create directory if doesn't exist
dir.create(loc_dir, showWarnings = FALSE)

# Save the data
fwrite(new_data, paste0(loc_dir, filename))

# Make a heatmap
if (do_heatplot) {
  ph1 <- pheatmap::pheatmap(new_data,
    filename = paste0(loc_dir, ph_save_name),
    main = "Heatmap of generated data"
  )
}

if (do_cor_heatplot) {
  ph2 <- pheatmap::pheatmap(cor(t(new_data)),
    filename = paste0(loc_dir, cor_ph_save_name),
    main = "Heatmap of correlation across generated data"
  )
}

# === Paul's stuff ============================================================

if (F) {
  univariateData <- c(
    rnorm(n_clust[1], means[1]),
    rnorm(n_clust[2], means[2]),
    rnorm(n_clust[3], means[3]),
    rnorm(n_clust[4], means[4]),
    rnorm(n_clust[5], means[5])
  )

  # ks <- density(univariateData)
  # plot(ks)


  # hist(univariateData, p)
  mvData <- c()
  for (i in 1:p) {
    mvData <- rbind(
      mvData,
      c(
        rnorm(n_clust[1], means[1]),
        rnorm(n_clust[2], means[2]),
        rnorm(n_clust[3], means[3]),
        rnorm(n_clust[4], means[4]),
        rnorm(n_clust[5], means[5])
      )
    )
  }

  pheatmap::pheatmap(t(mvData))
  pheatmap::pheatmap(t(mvData[1:100, ]))

  mvData2 <- c()
  for (i in 1:p) {
    mvData2 <- rbind(mvData2, univariateData + rnorm(n, mu_peturbation, sd_peturbation))
  }
}
