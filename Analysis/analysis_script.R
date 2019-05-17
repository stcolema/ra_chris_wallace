#!/usr/bin/env Rscript

# Example of a call:
# Rscript ra_chris_wallace/Data_prep/matrix_transposer.R -a TRUE
# -d ./MDI/Data/Original\ data/ -w ./VSN_NA_min_data/ -n TRUE
# --na_people_threshold 0.1 --na_probe_threshold 0.1 -t TRUE -v TRUE


# Rscript to convert data from CEDAR cohorts (found here: http://cedar-web.giga.ulg.ac.be/)
# Call from the command line with input arguments
# Writes to current directory

# === Libraries ================================================================

# MDI specific codes (not sure if still used)
source("/home/MINTS/sdc56/Desktop/MDI/mdipp-1.0.1/scripts/analysis.R") # install.packages("mcclust", dep = T)

# For posterioir similarity matrices
Rcpp::sourceCpp("/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/posterior_sim_mat.cpp") # install.packages("Rcpp", dep = T)

function_dir <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/Analysis_script_functions/"

function_scripts <- c(
  "plot_phi_series.R",
  "plot_phi_densities.R",
  "plot_phi_histograms.R",
  "plot_similarity_matrices.R",
  "plot_rand_index.R",
  "plot_fused_genes.R"
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

# For Bayesian clustering output (possibly should use mcclust) and Adjusted Rand Index
library(mclust) # install.packages("mclust", dep = T)

# For plotting
library(ggplot2)

# For command line arguments
library(optparse) # install.packages("optparse")

# Load data.table to access fread and fwrite functions
library(data.table) # install.packages("data.table", dep = T)

# Load magrittr for the pipe %>%
library(magrittr)

# For defensive programming in plot_fused_genes
library(attempt)

# === Functions ================================================================


# Function to find the mode of a vector
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Function to plot and save tree for a given similarity matrix
plot_tree <- function(sim_mat, plot_name, plot_type = ".png") {
  if (!(plot_type %in% c(".pdf", ".png"))) {
    stop("Wrong plot type attempted: must be one of '.pdf' or '.png'")
  }

  # Open a file
  if (plot_type == ".png") {
    png(plot_name)
  } else {
    pdf(plot_name)
  }
  # 2. Create a plot
  sim_mat %>%
    dist() %>%
    hclust() %>%
    plot()

  # Close the pdf file
  dev.off()
}

# For creating unique gene names
unique_names_from_recurring <- function(name_vec) {

  # Number of entries
  n <- length(name_vec)

  # Duplicate the name vec - this will hold the transformed versions
  duplicate_vec <- name_vec

  for (i in 1:n) {
    n_occurences <- sum(name_vec == name_vec[i])

    if (n_occurences > 1) {
      duplicate_vec[i] <- paste0(
        duplicate_vec[i],
        ".",
        sum(name_vec[1:i] == name_vec[i])
      )
    }
  }
  duplicate_vec
}


# User inputs from command line
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

    # Number of genes in each dataset
    optparse::make_option(c("--n_genes"),
      type = "integer",
      default = NA,
      help = "Number of genes in the datasets (if NA will find - this is included to avoid repetitive calculations) [default= %default]",
      metavar = "integer"
    ),

    # Number of iterations in MCMC step of MDI
    optparse::make_option(c("-n", "--n_iter"),
      type = "integer",
      default = NA,
      help = "Number of iterations in the MCMC step of MDI [default= %default]",
      metavar = "integer"
    ),

    # Thinning factor in MCMC step of MDI
    optparse::make_option(c("-t", "--thin"),
      type = "integer",
      default = 1,
      help = "Thinnging factor in the MCMC step of MDI [default= %default]",
      metavar = "integer"
    ),

    # Burn-in for MCMC step of MDI
    optparse::make_option(c("-b", "--burn"),
      type = "integer",
      default = 0,
      help = "Burn-in for the MCMC step of MDI [default= %default]",
      metavar = "integer"
    ),

    # .csv file describing which probes are not present in which datasets
    optparse::make_option(c("--probe_present"),
      type = "character",
      default = "/home/MINTS/sdc56/Desktop/MDI_small_geneset_outputs/Meta_data/probes_present_per_dataset.csv",
      help = "Name of .csv files describing which probes are empty in which dataset [default= %default]",
      metavar = "character"
    ),

    # .csv connecting probe IDs to genes
    optparse::make_option(c("--probe_key"),
      type = "character",
      default = "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/probe_key.csv",
      help = "Name of .csv files describing which probes are associated with which gene [default= %default]",
      metavar = "character"
    ),

    # Plot types to generate
    optparse::make_option(c("--plot_type"),
      type = "character",
      default = ".png",
      help = "Plot types, one of ``.png'' or ``.pdf'' [default= %default]",
      metavar = "character"
    ),

    # Instruction to plot trees
    optparse::make_option(c("--plot_trees"),
      type = "logical",
      default = FALSE,
      help = "Instruction to plot dendrogram from clustering on PSMs [default= %default]",
      metavar = "logical"
    ),

    # Instruction to plot adjusted rand index
    optparse::make_option(c("--plot_rand_index"),
      type = "logical",
      default = TRUE,
      help = "Instruction to plot adjusted rand index comparing clusterings [default= %default]",
      metavar = "logical"
    ),

    # Instruction to plot trees
    optparse::make_option(c("--plot_similarity_matrices"),
      type = "logical",
      default = TRUE,
      help = "Instruction to plot heatmaps of the clustering over iterations [default= %default]",
      metavar = "logical"
    ),

    # Instruction to plot phi scatter plots
    optparse::make_option(c("--plot_phi_series"),
      type = "logical",
      default = TRUE,
      help = "Instruction to plot phi parameters over MCMC iterations for each dataset [default= %default]",
      metavar = "logical"
    ),

    # Instruction to plot phi densities
    optparse::make_option(c("--plot_phi_densities"),
      type = "logical",
      default = TRUE,
      help = "Instruction to plot density of phi parameters for each dataset [default= %default]",
      metavar = "logical"
    ),

    # Instruction to plot phi densities
    optparse::make_option(c("--plot_phi_histograms"),
      type = "logical",
      default = TRUE,
      help = "Instruction to plot histogram of phi parameters for each dataset [default= %default]",
      metavar = "logical"
    ),

    # Location for expression data
    optparse::make_option(c("--expression_dir"),
      type = "character",
      default = "/home/MINTS/sdc56/Desktop/Matlab_input_small_names",
      help = "Path to the directory containing the expression data .csv files [default= %default]",
      metavar = "character"
    ),


    # Location for expression data
    optparse::make_option(c("--plot_expression_data"),
      type = "logical",
      default = TRUE,
      help = "Instruction to plot heatmaps of the expression data [default= %default]",
      metavar = "logical"
    ),

    # Instruction to time programme
    optparse::make_option(c("--seed"),
      type = "integer",
      default = 1,
      help = "Random seed for script [default= %default]",
      metavar = "integer"
    ),


    # Threshold for proportion of iterations we reuiqre probes to have the same
    # label to be considered ''fused''
    optparse::make_option(c("--fusion_threshold"),
      type = "double",
      default = 0.5,
      help = "Threshold for proportion of iterations we reuiqre probes to have the same label to be considered ''fused'' [default= %default]",
      metavar = "double"
    ),


    # Instruction to plot gene expression data hheatmaps for fused genes
    optparse::make_option(c("--plot_fused_genes"),
      type = "logical",
      default = TRUE,
      help = "Instruction to plot heatmaps of the expression data for fused genes [default= %default]",
      metavar = "logical"
    ),

    # Instruction to time programme
    optparse::make_option(c("--time"),
      type = "logical",
      default = FALSE,
      help = "instruciton to record runtime of function [default= %default]",
      metavar = "logical"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}


# === Main script ==============================================================

args <- input_arguments()
save_plots <- T
stm_i <- Sys.time()

seed <- args$seed
set.seed(seed)

# Read in the probes present - this is a matrix of bools with the first column
# as the probe IDs and the remaining columns corrresponding to the cell types
# with TRUE indicating the probe is present in this cell type (i.e. not added
# manually with an imputed value) and FALSE indicates we added it in.
probes_present_dt <- fread(args$probe_present)

# Read in the file relating the probe IDs to the related gene
probe_key <- fread(args$probe_key)

# Read in the MDI output file
file_path <- args$dir
# file_path <- "Analysis/MDI_runs/vsn_many_seeds/"

# MDI call specific values
n_iter <- args$n_iter
thin <- args$thin
burn <- args$burn

# Number of genes present in datasets
# The number of genes is the number of columns in the mcmc_output
# exclusing the columns containing information on the phi and mass parameters
n_genes <- args$n_genes

# The files / tissues used in the MDI
files_present <- args$datasets %>%
  strsplit(., " ") %>%
  unlist()

num_datasets <- length(files_present)
col_names <- paste0("D", 1:num_datasets)

# Remove the file extension
dataset_names <- tools::file_path_sans_ext(files_present)

do_dendrograms_ie_trees <- args$plot_trees
do_rand_plot <- args$plot_rand_index
do_similarity_matrices_plot <- args$plot_similarity_matrices
do_phis_series <- args$plot_phi_series
do_phis_densities <- args$plot_phi_densities
do_expression_heatmap <- args$plot_expression_data
do_phis_histograms <- args$plot_phi_histograms
do_fused_gene_expression <- args$plot_fused_genes

# Plto type
plot_type <- args$plot_type

# Fusion threshold
fusion_threshold <- args$fusion_threshold

# Create the common save path
save_path <- file_path

# For Rand index plots
# Save the common part of each title
generic_title <- "MDI: Adjusted Rand index for"
generic_save_name <- file_path


mdi_output_files <- list.files(path = file_path, full.names = T, include.dirs = F) %>%
  grep("csv", ., value = TRUE)

num_files <- length(mdi_output_files)
file_names <- basename(tools::file_path_sans_ext(mdi_output_files))

# Declare empty variable to capture information
mdi_allocation <- list()
allocation_list <- list()

eff_n_iter <- n_iter / thin # - burn

# For plotting phis
phis <- list()
count <- 0
start_index <- ceiling(burn / thin) + 1 # Consider letting this equal "burn"

# Convert this into useable output using functions provided by Sam Mason
mcmc_out_lst <- list() # vector("list", num_files)

# Large palette of qualitative colours
# (see https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vector <- brewer.pal(n = 12, name = "Paired") # unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

for (i in 1:num_files) {
  curr_name <- file_names[[i]]
  # comp_mcmc_out <- readMdiMcmcOutput(mdi_output_files_2)
  # mcmc_out_lst[[curr_name]] <- readMdiMcmcOutput(mdi_output_files[[i]])
  # my_mcmc_out <-  readMdiMcmcOutput(mdi_output_files[[i]])
  mcmc_out_lst[[curr_name]] <- fread(mdi_output_files[[i]])
}


# mcmc_output <- readMdiMcmcOutput(mdi_output_file)

# clust_occ <- mcmc_out_lst %>%
# lapply(getClustersOccupied)

# clust_occ <- getClustersOccupied(mcmc_output)
# head(clust_occ)

# By the 8 x 25 = 200 iteration we have reached stationarity in the number of
# clusters occupied
# summary(clust_occ[8:800, ])

# NOte that as each mdi output is on the same data, we only need to count one dataset
if (is.na(n_genes)) {
  n_genes <- mcmc_out_lst[[1]] %>%
    select(contains("Dataset1")) %>%
    ncol()

  # n_genes <- mcmc_out_lst[[1]]$nitems
}

# === Plotting phis ==========================================================

if (do_phis_series) {
  plot_phi_series(mcmc_out_lst,
    file_path,
    num_files,
    num_datasets,
    start_index,
    eff_n_iter,
    save_plots = T
  )

  # print("Saving plots of phi as a time series.")
  #
  # plot_phi_series(
  #   mcmc_out_lst,
  #   file_path,
  #   num_files,
  #   num_datasets,
  #   start_index,
  #   eff_n_iter,
  #   save_plots
  # )
  #
  # # Create directory to save this output in
  # loc_dir <- paste0(file_path, "Phi_series_plots/")
  # dir.create(loc_dir, showWarnings = FALSE)
  #
  # # Iterate over the files and then the combinations of phis
  # for (i in 1:num_files) {
  #
  #   # For heatmapping the phis across datasets
  #   phi_comparison_df <- as.data.frame(matrix(
  #     nrow = num_datasets,
  #     ncol = num_datasets,
  #     0
  #   )) %>%
  #     magrittr::set_colnames(dataset_names) %>%
  #     magrittr::set_rownames(dataset_names)
  #
  #   for (j in 1:(num_datasets - 1)) {
  #     for (k in (j + 1):num_datasets) {
  #
  #       # Count for the index of the list object where phis are stored
  #       count <- count + 1
  #
  #       col_name <- paste0("Phi_", j, k)
  #
  #       # Pull out the column for the relevant phi
  #       phis[[count]] <- mcmc_out_lst[[i]] %>%
  #         select(contains(col_name))
  #
  #       # Which tissues
  #       dataset_j <- dataset_names[[j]]
  #       dataset_k <- dataset_names[[k]]
  #
  #       # Create plot labels
  #       plot_title <- bquote(Phi ~ "for" ~ .(dataset_j) ~ "and" ~ .(dataset_k))
  #       y_axis_title <- substitute(Phi[ind1], list(ind1 = paste0(j, k)))
  #       sub_title <- paste("Iterations", (burn + thin), "through", n_iter)
  #
  #       # The save file name
  #       save_title <- paste0(loc_dir, "file_", i, "_Phi_", j, k, plot_type)
  #
  #       if (save_plots) {
  #         # Open graphic to save plot to
  #         if (plot_type == ".pdf") {
  #           pdf(save_title)
  #         } else {
  #           png(save_title)
  #         }
  #       }
  #
  #       # cat("Start index:" )
  #       # print(start_index)
  #       # cat("Number of iterations:")
  #       # print(n_iter)
  #       #
  #       # cat("Effective number of iterations:")
  #       # print(eff_n_iter)
  #       #
  #       # print("Description of phi parameter")
  #       # print(str(phis[[count]][[1]]))
  #       #
  #       # print("That's phi")
  #       # print("")
  #
  #       # Plot Phi vs index (ignore initial values as misleading)
  #       plot(start_index:eff_n_iter, phis[[count]][[1]][start_index:eff_n_iter],
  #         main = plot_title,
  #         col.main = "black",
  #         sub = sub_title,
  #         col.sub = "blue",
  #         ylab = y_axis_title,
  #         xlab = "Index"
  #       )
  #
  #       # Close graphic
  #       if (save_plots) {
  #         dev.off()
  #       }
  #
  #       mean_phi <- phis[[count]][[1]][start_index:eff_n_iter] %>% mean()
  #
  #       phi_comparison_df[j, k] <- phi_comparison_df[k, j] <- mean_phi
  #     }
  #   }
  #
  #   # Heatmap the average phi (after some burn in)
  #   phi_pheatmap_title <- "Heatmap comparing phis across datasets"
  #   phi_pheatmap_file_name <- paste0(save_path, "Phi_heatmap_", i, plot_type)
  #
  #   if (save_plots) {
  #     pheatmap(phi_comparison_df,
  #       main = phi_pheatmap_title,
  #       filename = phi_pheatmap_file_name
  #     )
  #
  #     # dev.off()
  #   } else {
  #     pheatmap(phi_comparison_df,
  #       main = phi_pheatmap_title,
  #       filename = phi_pheatmap_file_name
  #     )
  #   }
  # }
}

# === Prepare the tibble =======================================================

# Create an empty dataframe with column names corresponding to dataset numbers
compare_df <- data.frame(matrix(ncol = num_datasets, nrow = n_genes))
# colnames(compare_df) <- col_names

# Remove the output file from this list if you saved it in the same location as
# the data files
# if (mdi_output_file %in% files_present) {
#   index_to_empty <- files_present == mdi_output_file
#   files_present %<>%
#     list.remove(index_to_empty)
# }

# Have to have the order from the call
# CD4, CD8, CVD19, CD14, CD15, platelets, ileonic, colonic and rectal biopsies
# transcriptome data for six circulating immune cell types (CD4+ T lymphocytes,
# CD8+ T lymphocytes, CD19+ B lymphocytes, CD14+ monocytes, CD15+ granulocytes,
# platelets) as well as ileal, colonic, and rectal biopsies (IL, TR, RE)
# 323 healthy Europeans
# files_present <- c("CD14.csv", "CD15.csv", "IL.csv", "CD4.csv", "CD8.csv", "CD19.csv", "PLA.csv", "RE.csv", "TR.csv")
# dataset_names <- tools::file_path_sans_ext(files_present)

# Set the cell type to the column names (make sure order is as per MDI call)
colnames(compare_df) <- dataset_names

print("Constructing tibble.")

# Now put everything in a tibble
compare_tibble <- tibble(
  mdi = rep(file_names, num_datasets),
  dataset = rep(dataset_names, num_files), # unlist(lapply(dataset_names, rep, num_files))
  seed = unlist(lapply(1:num_files, rep, num_datasets)), # rep(1:num_files, num_datasets),
  pred_allocation = list(compare_df),
  mdi_allocation = rep(vector("list", num_files), num_datasets),
  similarity_matrix = list(data.frame(matrix(ncol = n_genes, nrow = n_genes))),
  expression_data = rep(vector("list", num_files), num_datasets),
  non_zero_probes_ind = rep(vector("list", num_files), num_datasets),
  non_zero_probes = rep(vector("list", num_files), num_datasets),
  fused_probes = rep(vector("list", num_files), num_datasets),
  fused_non_zero_probes = rep(vector("list", num_files), num_datasets)
)

# === MDI output ===============================================================

# Global similarity between datasets
phis <- list()

# x4 <- sample(1:eff_iter, 50, replace=F)

for (j in 1:num_files) {
  # compare_tibble$MDI_allocation[j] <- list()
  # compare_tibble$Pred_allocation[j] <- list()

  mdi_pos_sim_mat <- list()
  check_median_makes_sense_map <- list()

  # Capture the allocation information in the named lists and the predicted
  # allocation in the dataframe

  phis[[j]] <- mcmc_out_lst[[j]] %>%
    dplyr::select(contains("Phi"))


  for (i in 1:num_datasets) {
    dataset_name <- paste0("Dataset", i)

    # Get the allocation and drop the burn in
    # mdi_allocation[[i]] <- .mdi_alloc <- getMdiAllocations(mcmc_out_lst[[j]], i)

    # mdi_alloc2 <- getMdiAllocations(comp_mcmc_out, i)

    mdi_allocation[[i]] <- .mdi_alloc <- mcmc_out_lst[[j]] %>%
      dplyr::select(contains(dataset_name)) %>%
      magrittr::extract(start_index:eff_n_iter, )

    # mdi_pos_sim_mat[[i]][[1]] <- .sim_mat <- similarity_mat(t(.mdi_alloc[,1:1000]))

    # print(str(.mdi_alloc))
    # print(start_index)
    # print(eff_n_iter)

    mdi_pos_sim_mat[[i]] <- .sim_mat <- similarity_mat(t(.mdi_alloc))

    # print("This point.")

    if (i == 1) {
      probe_names <- colnames(.mdi_alloc) %>%
        stringr::str_remove_all(paste0(dataset_name, "_")) %>%
        stringr::str_remove_all("X")
    }

    row.names(.sim_mat) <- colnames(.sim_mat) <- probe_names

    # print("Now over here")

    colnames(.mdi_alloc) <- probe_names

    compare_tibble$mdi_allocation[i + (j - 1) * num_files][[1]] <- .mdi_alloc

    compare_tibble$similarity_matrix[i + (j - 1) * num_files][[1]] <- .sim_mat

    # sams_pos_sim <- genPosteriourSimilarityMatrix(.mdi_alloc)
    #
    # # Individuals are columns rather than rows (henc transpose)
    # mdi_pos_sim_mat[[i]] <- list()
    # mdi_pos_sim_mat[[i]][[1]] <- .sim_mat <- similarity_mat(t(.mdi_alloc[,1:1000]))
    # mdi_pos_sim_mat[[i]][[2]] <- .sim_mat <- similarity_mat(t(.mdi_alloc[,5000:6000]))
    # mdi_pos_sim_mat[[i]][[3]] <- .sim_mat <- similarity_mat(t(.mdi_alloc[,11000:12000 ]))
    # # break
    # # By checking the imilarity of each row we can decide if the median is an
    # # accurate method to allocate class (if all rows are highly similar label
    # # flipping did not occur)
    # check_median_makes_sense_map[[i]] <- .sense_check <- similarity_mat(.mdi_alloc[c(seq(1, eff_n_iter, 25), eff_n_iter), ])
    # break
    # }

    # pheatmap(mdi_pos_sim_mat[[1]])
    # pheatmap(.sim_mat)
    # x1 <- hclust(dist(.sim_mat), method = "complete", members = NULL)
    # plot(x1)


    allocation_list[[i]] <- .pred_alloc <- .sim_mat %>%
      mcclust::Simtocl()
    # dist() %>%
    # hclust() %>%
    # cutree(h = 3.2) # abritrary number here based on pheatmap of .sim_mat

    # allocation_list[[i]] <- .pred_alloc <- apply(.mdi_alloc[-(1:burn), ], 2, getmode)
    # compare_df[[dataset_names[i]]] <- .pred_alloc

    # print(.pred_alloc)

    compare_tibble$pred_allocation[i + (j - 1) * num_files][[1]] <- .pred_alloc


    # if (i == 1) {
    # row.names(compare_df) <- names(.pred_alloc)
    # }
  }
  # compare_tibble$MDI_allocation[j] <- mdi_allocation
  # compare_tibble$Pred_allocation[j] <- compare_df
}

# saveRDS(
#   compare_tibble,
#   paste0(save_path, "compare_tibble.rds")
# )

# === Find GENE names ==========================================================

probe_key_rel <- probe_key[probe_key$ProbeID %in% probe_names, ]

# Pull out the Gene IDs in the correct order
gene_id <- probe_key_rel %>%
  .[match(probe_names, .$ProbeID)] %>%
  .$Gene

unique_gene_id <- gene_id %>%
  unique_names_from_recurring()

# === Plot PSM trees ===========================================================
if (do_dendrograms_ie_trees) {
  print("Saving dendrogram plots.")

  loc_dir <- paste0(file_path, "Tree_plots/")
  dir.create(loc_dir, showWarnings = FALSE)

  # Plot the dendrograms for the PSM
  for (i in 1:nrow(compare_tibble)) {
    dataset <- dataset_names[[i]]
    plot_name <- paste0(save_path, dataset, "_tree", plot_type)
    plot_tree(compare_tibble$similarity_matrix[i][[1]], plot_name)
  }
}


# === Probes present ===========================================================

# Find which probes are relevant from the full set
probes_actually_present_ind <- probes_present_dt %>%
  magrittr::use_series("V1") %>%
  magrittr::is_in(probe_names)

# Select these
probes_actually_present <- probes_present_dt[probes_actually_present_ind, ]

# Find the appropriate order
probes_order <- match(probe_names, probes_actually_present$V1)

# Order the probes so comparable to allocation data frame
probes_present_ordered <- probes_actually_present[probes_order, ] %>%
  as.data.frame()

# Remove the irrelevant columns and set row names
probes_present_final <- probes_present_ordered %>%
  magrittr::set_rownames(probes_actually_present$V1) %>%
  magrittr::extract(, -c(1, 3, 8))

colnames(probes_present_final) <- dataset_names

# === Phi denisty plots ================================================================

if (do_phis_densities) {
  print("Saving phi density plots.")
  plot_phi_densities(phis, file_path)

  # loc_dir <- paste0(save_path, "Phi_density_plots/")
  # dir.create(loc_dir, showWarnings = FALSE)
  #
  # for (i in 1:length(phis)) {
  #   for (j in 1:ncol(phis[[i]])) {
  #     curr_phi <- colnames(phis[[i]])[j]
  #     # plot_name <- paste0(loc_dir, curr_phi, "_plot", plot_type)
  #
  #     # png(plot_name)
  #     #
  #     # phis[[i]][, ..j] %>%
  #     #   unlist() %>%
  #     #   plot()
  #     #
  #     # dev.off()
  #
  #     density_plot_name <- paste0(loc_dir, curr_phi, "_density_plot", plot_type)
  #
  #     density_title <- paste(curr_phi, ": density plot")
  #
  #     # Find which datasets are used here
  #     dataset_indices <- readr::parse_number(curr_phi)
  #
  #     # The following steps do not work if we have more than 9 datasets (which is not relevant to me)
  #     # This extracts the first dataset index (i.e. from 67 takes 6)
  #     dataset_index_1 <- dataset_indices[1] %>%
  #       "/"(10) %>%
  #       floor(.)
  #
  #     # This extracts the secodn dataset index
  #     dataset_index_2 <- dataset_indices[1] %>% "%%"(10)
  #
  #     dataset_1 <- dataset_names[dataset_index_1]
  #     dataset_2 <- dataset_names[dataset_index_2]
  #
  #     density_subtitle <- paste(
  #       "Iterations",
  #       (burn + thin),
  #       "through",
  #       n_iter,
  #       "for",
  #       dataset_1,
  #       "and",
  #       dataset_2
  #     )
  #
  #     ggplot(data = phis[[1]][start_index:n_iter, ], aes_string(x = curr_phi)) +
  #       geom_density() +
  #       labs(
  #         title = density_title,
  #         subtitle = density_subtitle
  #       )
  #
  #     ggsave(density_plot_name)
  #   }
  # }
}

# === Phi histograms ===========================================================

if (do_phis_histograms) {
  print("Saving phi histogram plots.")
  plot_phi_histograms(phis, file_path)

  # loc_dir <- paste0(save_path, "Phi_histograms/")
  # dir.create(loc_dir, showWarnings = FALSE)
  #
  # for (i in 1:length(phis)) {
  #   for (j in 1:ncol(phis[[i]])) {
  #     curr_phi <- colnames(phis[[i]])[j]
  #     plot_name <- paste0(loc_dir, curr_phi, "_plot", plot_type)
  #
  #
  #     # Find which datasets are used here
  #     dataset_indices <- readr::parse_number(curr_phi)
  #
  #     # The following steps do not work if we have more than 9 datasets (which is not relevant to me)
  #     # This extracts the first dataset index (i.e. from 67 takes 6)
  #     dataset_index_1 <- dataset_indices[1] %>%
  #       "/"(10) %>%
  #       floor(.)
  #
  #     # This extracts the secodn dataset index
  #     dataset_index_2 <- dataset_indices[1] %>% "%%"(10)
  #
  #     dataset_1 <- dataset_names[dataset_index_1]
  #     dataset_2 <- dataset_names[dataset_index_2]
  #
  #     histogram_plot_name <- paste0(loc_dir, curr_phi, "_histogram_plot", plot_type)
  #
  #     histogram_title <- paste(curr_phi, ": histogram plot")
  #
  #     histogram_subtitle <- paste(
  #       "Iterations",
  #       (burn + thin),
  #       "through",
  #       n_iter,
  #       "for",
  #       dataset_1,
  #       "and",
  #       dataset_2
  #     )
  #
  #     ggplot(data = phis[[1]][start_index:n_iter, ], aes_string(x = curr_phi)) +
  #       geom_histogram() +
  #       labs(
  #         title = histogram_title,
  #         subtitle = histogram_subtitle
  #       )
  #
  #     ggsave(histogram_plot_name)
  #   }
  # }
}


# ggplot(data = phis[[1]], aes(x = Phi_12)) +
#   geom_density()


# === Plot posterior similarity matrices =======================================

# If making heatplots of the clusterings across iterations
if (do_similarity_matrices_plot) {
  print("Saving heatmaps of PSMs.")

  plot_similarity_matrices(
    compare_tibble$similarity_matrix,
    probes_present_final,
    dataset_names,
    num_files,
    num_datasets,
    file_path
  )


  # loc_dir <- paste0(save_path, "Similarity_matrices/")
  # dir.create(loc_dir, showWarnings = FALSE)
  #
  # # Define the type of file to save
  # # plot_type <- ".pdf" # ".png" or ".pdf"
  #
  # for (j in 1:num_files) {
  #   # curr_save_path <- paste0(save_path, "seed_", j)
  #   # Iterate over datasets
  #   for (i in 1:num_datasets) {
  #     # Find the dataset name
  #     dataset <- dataset_names[[i]]
  #
  #     # Create the save location and file name
  #     file_name <- paste0(loc_dir, "similarity_matrix_", dataset, plot_type)
  #
  #     # Create the title of the plot
  #     title <- paste("Similarity matrix for", dataset)
  #
  #     rel_rows <- probes_present_final[, i]
  #
  #     # re_mat <- .sim_mat
  #     # re_mat[! rel_rows, ! rel_rows] <- 100
  #     #
  #     # ph1 <- pheatmap(re_mat)
  #     # ord1 <- ph1$tree_col[["order"]]
  #     # pheatmap(.sim_mat[ord1, ord1], cluster_cols = F, cluster_rows = F)
  #     #
  #     # rel_sim_mat <- compare_tibble$similarity_matrix[i + (j - 1) * num_files][[1]] %>%
  #     #   magrittr::extract(., rel_rows, )
  #
  #     # Make the heatmap (note that we do not cluster rows).
  #     # We use a sample of the iterations as otherwise it becomes too heavy
  #
  #     # print(mcclust::Simtocl(compare_tibble$similarity_matrix[i + (j - 1) * num_files][[1]][1:4, 1:4]))
  #     # print(compare_tibble$similarity_matrix[i + (j - 1) * num_files][[1]][1:4, 1:4])
  #
  #     if (save_plots) {
  #       ph <- pheatmap(compare_tibble$similarity_matrix[i + (j - 1) * num_files][[1]],
  #         main = title,
  #         cluster_rows = T,
  #         filename = file_name
  #       )
  #     } else {
  #       ph <- pheatmap(compare_tibble$similarity_matrix[i + (j - 1) * num_files][[1]],
  #         main = title,
  #         cluster_rows = T
  #       )
  #     }
  #   }
  # }
}
# pheatmap::pheatmap(check_median_makes_sense_map[[1]], cluster_rows = F, cluster_cols = T)
# pheatmap::pheatmap(check_median_makes_sense_map[[2]], cluster_rows = F, cluster_cols = T)
# pheatmap::pheatmap(check_median_makes_sense_map[[3]], cluster_rows = F, cluster_cols = T)
# pheatmap::pheatmap(check_median_makes_sense_map[[4]], cluster_rows = F, cluster_cols = T)
# pheatmap::pheatmap(check_median_makes_sense_map[[5]], cluster_rows = F, cluster_cols = T)
# pheatmap::pheatmap(check_median_makes_sense_map[[6]], cluster_rows = F, cluster_cols = T)
# pheatmap::pheatmap(check_median_makes_sense_map[[7]], cluster_rows = F, cluster_cols = T)
# pheatmap::pheatmap(check_median_makes_sense_map[[8]], cluster_rows = F, cluster_cols = T)
# pheatmap::pheatmap(check_median_makes_sense_map[[9]], cluster_rows = F, cluster_cols = T)

# === Plot adjusted rand index==================================================
# If instructed to make Rand index plots
if (do_rand_plot) {
  print("Saving scatter plots of adjusted rand index comparing final clustering to clustering at each iteration.")

  plot_rand_index(
    compare_tibble$mdi_allocation,
    file_path,
    dataset_names,
    num_files,
    num_datasets,
    eff_n_iter,
    burn
  )


  # loc_dir <- paste0(save_path, "Adjusted_rand_index_plots/")
  # dir.create(loc_dir, showWarnings = FALSE)
  #
  # # Define the plot type
  # # plot_type <- ".pdf" # ".png" or ".pdf"
  #
  # # Create the lists to hold the dataset specific results
  # rand <- list()
  # rand_plots <- list()
  #
  #
  # for (j in 1:num_files) {
  #
  #   # Loop over the datasets
  #   for (i in 1:num_datasets) {
  #     # The current dataset name
  #     dataset <- dataset_names[[i]]
  #
  #     # Append this to the generic title to create the specific title
  #     curr_title <- paste(generic_title, dataset)
  #
  #     # Similarly for the name of the save file
  #     curr_save_file <- paste0(loc_dir, "rand_index_plot_", dataset, plot_type)
  #
  #     # Create a vector of the adjusted rand index comparing the modal cluster
  #     # to the clustering at each iteration
  #     rand[[i + (j - 1) * num_files]] <- apply(
  #       compare_tibble$mdi_allocation[i + (j - 1) * num_files][[1]],
  #       # mdi_allocation[[i]],
  #       1,
  #       # mcclust::arandi,
  #       mclust::adjustedRandIndex,
  #       compare_tibble$mdi_allocation[i][[1]][eff_n_iter - burn, ]
  #       # compare_df[, i]
  #     )
  #
  #     # As we use ggplot2, put this in a dataframe
  #     plot_data <- data.frame(Rand_index = rand[[i]], Iteration = 1:length(rand[[i]]))
  #
  #     # Plot the Rand index against iteration number
  #     rand_plots[[i]] <- ggplot2::ggplot(data = plot_data) +
  #       ggplot2::geom_point(ggplot2::aes(x = Iteration, y = Rand_index)) +
  #       ggplot2::labs(
  #         title = curr_title,
  #         subtitle = "Comparing clustering in last iteration to clustering at each iteration",
  #         x = "Index",
  #         y = "Adjusted Rand Index"
  #       )
  #
  #     if (save_plots) {
  #       # Save the plot
  #       ggplot2::ggsave(curr_save_file, plot = rand_plots[[i]])
  #     }
  #   }
  # }
}

# === Heatmap expression data ==================================================
loc_dir <- paste0(save_path, "Expression_heatmaps/")

if (do_expression_heatmap) {
  print("Saving gene expression heatmaps.")


  dir.create(loc_dir, showWarnings = FALSE)
}


# Directory holding the expression data files
data_dir <- args$expression_dir

# Generic title and filename for pheatmap
gen_ph_title <- ": heatmap of expression data"
gen_ph_file_name <- paste0(loc_dir, "pheatmap_")

# Find the expression data files
mdi_input_files <- list.files(path = data_dir, full.names = T, include.dirs = F) %>%
  grep("csv", ., value = TRUE)

input_file_names <- basename(tools::file_path_sans_ext(mdi_input_files))

expression_datasets <- input_file_names %>%
  gsub("([^\\_]+)\\_.*", "\\1", .)
# stringr::str_replace("_sma_mat_nv", "")

datasets_relevant_indices <- files_present %>%
  tools::file_path_sans_ext() %>%
  match(expression_datasets)

datasets_relevant <- expression_datasets[datasets_relevant_indices]
relevant_input_files <- mdi_input_files[datasets_relevant_indices]

mega_df <- data.frame(matrix(nrow = n_genes, ncol = 0)) %>%
  magrittr::set_rownames(probe_names)

big_annotation <- data.frame(matrix(nrow = n_genes, ncol = num_datasets)) %>%
  magrittr::set_rownames(probe_names) %>%
  magrittr::set_colnames(datasets_relevant)

n_total_clusters <- 0

data_files <- list()
for (i in 1:num_datasets) {
  curr_dataset <- datasets_relevant[[i]]

  file_name <- gen_ph_file_name %>%
    paste0(datasets_relevant[[i]], plot_type)

  ph_title <- datasets_relevant[[i]] %>%
    paste0(gen_ph_title)

  # Prepare the clustering information as an annotation row for pheatmap
  pred_clustering <- compare_tibble$pred_allocation[compare_tibble$dataset == curr_dataset][[1]] %>%
    as.factor() %>%
    as.data.frame() %>%
    magrittr::set_rownames(probe_names) %>%
    magrittr::set_colnames(c("Cluster"))


  # Read in the expression data
  f <- relevant_input_files[[i]]
  expression_data <- fread(f)

  # Tidy (remove NAs and row name column) and convert to the appropriate format
  # for pheatmap
  expression_data_tidy <- expression_data %>%
    magrittr::extract(, -1) %>%
    as.matrix() %>%
    magrittr::set_rownames(probe_names)

  expression_data_tidy[is.na(expression_data_tidy)] <- 0

  # print(str(compare_tibble$expression_data[compare_tibble$dataset == curr_dataset]))

  # print(compare_tibble$non_zero_probes_ind) <- .non_zero_probes <- rowSums(expression_data_tidy) != 0)

  # Add the expression data to our tibble
  num_occurences_dataset <- length(compare_tibble$expression_data[compare_tibble$dataset == curr_dataset])
  for (k in 1:num_occurences_dataset) {
    compare_tibble$expression_data[compare_tibble$dataset == curr_dataset][[k]] <- expression_data_tidy
    compare_tibble$non_zero_probes_ind[compare_tibble$dataset == curr_dataset][[k]] <- .non_zero_probes <- rowSums(expression_data_tidy) != 0
    compare_tibble$non_zero_probes[compare_tibble$dataset == curr_dataset][[k]] <- names(.non_zero_probes)[.non_zero_probes]
  }

  data_files[[i]] <- expression_data

  # Specify colors based on cluster labels
  cluster_labels <- levels(pred_clustering$Cluster)
  n_clusters <- length(cluster_labels)

  n_total_clusters <- max(n_clusters, n_total_clusters)

  if (do_expression_heatmap) {
    if (n_clusters > 12) {
      print("Too many clusters. Cannot include annotation row.")

      # Pheatmap
      expression_data_tidy %>%
        pheatmap(
          filename = file_name,
          main = ph_title
        )
    } else {
      col_pal <- col_vector[1:n_clusters] %>%
        # sample(col_vector, n_clusters) %>%
        magrittr::set_names(cluster_labels)

      annotation_colors <- list(Cluster = col_pal)

      row_order <- pred_clustering$Cluster %>% order()

      # Pheatmap
      expression_data_tidy[row_order, ] %>%
        pheatmap(
          filename = file_name,
          main = ph_title,
          cluster_rows = F,
          annotation_row = pred_clustering,
          annotation_colors = annotation_colors
        )
    }
  }

  # Record the clustering for the current dataset for annotation purposes
  big_annotation[[curr_dataset]] <- as.factor(pred_clustering$Cluster)

  # Record the expression data from the current dataset in the compound dataset
  mega_df <- mega_df %>%
    dplyr::bind_cols(as.data.frame(expression_data_tidy))
}

# row.names(mega_df) <- probe_names

big_file_name <- gen_ph_file_name %>%
  paste0("all_datasets", plot_type)

big_ph_title <- "All datasets" %>%
  paste0(gen_ph_title)

mega_matrix <- mega_df %>%
  as.matrix() %>%
  magrittr::set_rownames(probe_names)

if (do_expression_heatmap) {
  if (TRUE) { # n_total_clusters > 20) {
    pheatmap(mega_matrix,
      filename = big_file_name,
      main = big_ph_title
    )
  } else {
    col_pal <- sample(col_vector, n_total_clusters) %>%
      magrittr::set_names(cluster_labels)

    annotation_colors <- list(Cluster = col_pal)

    pheatmap(mega_matrix,
      filename = big_file_name,
      main = big_ph_title,
      annotation_row = big_annotation,
      annotation_colors = annotation_colors
    )
  }
}

# === Fused probes =============================================================

# Find which probes are ''fused'' across datasets
# We save this as a named list to the tibble. Each entry in the list corresponds
# to a dataset and records the fused probes between the current dataset and the
# entry name
print("Finding ''fused'' probes.")
for (k in 1:num_files) {
  for (i in 1:num_datasets) {
    dataset_i <- dataset_names[i]

    fused_probes_curr <- vector("list", num_datasets) %>%
      magrittr::set_names(dataset_names)

    fused_non_zero_probes_curr <- vector("list", num_datasets) %>%
      magrittr::set_names(dataset_names)

    # fused_probes_curr[dataset_i] <- 1.0

    for (j in 1:num_datasets) {
      dataset_j <- dataset_names[j]

      count <- count + 1

      # Find the unnormalised count for the numebr of times each probe has the same label across iterations
      .fusion_count <- compare_tibble$mdi_allocation[compare_tibble$dataset == dataset_i][[k]] %>%
        magrittr::equals(compare_tibble$mdi_allocation[compare_tibble$dataset == dataset_j][[k]])

      # convert this to a propotion
      .fusion_prob <- (1 / nrow(.fusion_count)) * colSums(.fusion_count)

      # Record as ''fused'' those probes for which the proportion of times they
      # have a common labelling exceeds some user-defined threshold (defualt of 0.5)

      # compare_tibble$fused_probes[
      #   compare_tibble$dataset == dataset_i
      # ][[k]] <-

      # This is rather ugly - we are recording a named list associating each
      # dataset with each other one
      fused_probes_curr[[dataset_j]] <- .fused_probes_ind_i <- .fusion_prob > fusion_threshold

      # compare_tibble$fused_probes[
      #   compare_tibble$dataset == dataset_j
      #   ][[k]] <- .fused_probes_ind_j <- .fusion_prob > fusion_threshold

      # Now record the "fused" probes excluding those which have a value of 0 across the dataset
      # compare_tibble$fused_non_zero_probes[
      #   compare_tibble$dataset == dataset_i
      # ][[k]]

      fused_non_zero_probes_curr[[dataset_j]] <- .fused_probes_ind_i & compare_tibble$non_zero_probes_ind[compare_tibble$dataset == dataset_i][[k]]

      # compare_tibble$fused_non_zero_probes[
      #   compare_tibble$dataset == dataset_j
      # ][[k]] <- .fused_probes_ind_j & compare_tibble$non_zero_probes_ind[compare_tibble$dataset == dataset_j][[k]]
    }

    # Record the named lists within the tibble
    compare_tibble$fused_probes[
      compare_tibble$dataset == dataset_i
    ][[k]] <- fused_probes_curr

    compare_tibble$fused_non_zero_probes[
      compare_tibble$dataset == dataset_i
    ][[k]] <- fused_non_zero_probes_curr
  }
}

# === Gene expression for fused and unfused genes ==============================

if (do_fused_gene_expression) {
  
  print("Saving heatmaps of pairwise fused gene expression across datasets.")
  
  fused_gene_heatmaps(
    compare_tibble$expression_data,
    compare_tibble$fused_probes,
    unique_gene_id,
    dataset_names,
    file_path,
    num_datasets,
    plot_type
  )
}

# === Timing ===================================================================

saveRDS(
  compare_tibble,
  paste0(save_path, "compare_tibble.rds")
)

if (args$time) {
  print((Sys.time() - stm_i))
}
