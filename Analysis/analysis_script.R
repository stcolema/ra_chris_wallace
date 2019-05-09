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
    optparse::make_option(c("--plot_heatmap_clusterings"),
      type = "logical",
      default = TRUE,
      help = "Instruction to plot adjusted rand index comparing clusterings [default= %default]",
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

# Read in the probes present - this is a matrix of bools with the first column
# as the probe IDs and the remaining columns corrresponding to the cell types
# with TRUE indicating the probe is present in this cell type (i.e. not added
# manually with an imputed value) and FALSE indicates we added it in.
probes_present_dt <- fread(args$probe_present)

# Read in the file relating the probe IDs to  the related gene
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
do_heatplot_clusterings <- args$plot_heatmap_clusterings

# Plto type
plot_type <- args$plot_type

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
start_index <- burn + 1 # Consider letting this equal "burn"

# Convert this into useable output using functions provided by Sam Mason
mcmc_out_lst <- list() # vector("list", num_files)

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

# Iterate over the files and then the combinations of phis
for (i in 1:num_files) {

  # For heatmapping the phis across datasets
  phi_comparison_df <- as.data.frame(matrix(
    nrow = num_datasets,
    ncol = num_datasets,
    0
  )) %>%
    magrittr::set_colnames(dataset_names) %>%
    magrittr::set_rownames(dataset_names)

  for (j in 1:(num_datasets - 1)) {
    for (k in (j + 1):num_datasets) {

      # Count for the index of the list object where phis are stored
      count <- count + 1

      col_name <- paste0("Phi_", j, k)

      # Pull out the column for the relevant phi
      phis[[count]] <- mcmc_out_lst[[i]] %>%
        select(contains(col_name))

      # Which tissues
      dataset_j <- dataset_names[[j]]
      dataset_k <- dataset_names[[k]]

      # Create plot labels
      plot_title <- bquote(Phi ~ "for" ~ .(dataset_j) ~ "and" ~ .(dataset_k))
      y_axis_title <- substitute(Phi[ind1], list(ind1 = paste0(j, k)))
      sub_title <- paste("Iterations", start_index, "through", n_iter)

      # The save file name
      save_title <- paste0(file_path, "file_", i, "_Phi_", j, k, plot_type)

      if (save_plots) {
        # Open graphic to save plot to
        if (plot_type == ".pdf") {
          pdf(save_title)
        } else {
          png(save_title)
        }
      }

      # Plot Phi vs index (ignore initial values as misleading)
      plot(start_index:n_iter, phis[[count]][[1]][start_index:n_iter],
        main = plot_title,
        col.main = "black",
        sub = sub_title,
        col.sub = "blue",
        ylab = y_axis_title,
        xlab = "Index"
      )

      # Close graphic
      if (save_plots) {
        dev.off()
      }

      mean_phi <- phis[[count]][[1]][start_index:n_iter] %>% mean()

      phi_comparison_df[j, k] <- phi_comparison_df[k, j] <- mean_phi
    }
  }

  # Heatmap the average phi (after some burn in)
  phi_pheatmap_title <- "Heatmap comparing phis across datasets"
  phi_pheatmap_file_name <- paste0(save_path, "Phi_heatmap_", i, plot_type)

  if (save_plots) {
    pheatmap(phi_comparison_df,
      main = phi_pheatmap_title,
      filename = phi_pheatmap_file_name
    )

    # dev.off()
  } else {
    pheatmap(phi_comparison_df,
      main = phi_pheatmap_title,
      filename = phi_pheatmap_file_name
    )
  }
}

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

# Now put everything in a tibble
compare_tibble <- tibble(
  mdi = rep(file_names, num_datasets),
  dataset = rep(dataset_names, num_files), # unlist(lapply(dataset_names, rep, num_files))
  seed = unlist(lapply(1:num_files, rep, num_datasets)), # rep(1:num_files, num_datasets),
  pred_allocation = list(compare_df),
  mdi_allocation = rep(vector("list", num_files), num_datasets),
  similarity_matrix = list(data.frame(matrix(ncol = n_genes, nrow = n_genes)))
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
      dplyr::select(contains(dataset_name))

    # mdi_pos_sim_mat[[i]][[1]] <- .sim_mat <- similarity_mat(t(.mdi_alloc[,1:1000]))
    mdi_pos_sim_mat[[i]] <- .sim_mat <- similarity_mat(t(.mdi_alloc))


    if (i == 1) {
      probe_names <- colnames(.mdi_alloc) %>%
        stringr::str_remove_all(paste0(dataset_name, "_")) %>%
        stringr::str_remove_all("X")
    }

    row.names(.sim_mat) <- colnames(.sim_mat) <- probe_names

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
      dist() %>%
      hclust() %>%
      cutree(h = 3.2) # abritrary number here based on pheatmap of .sim_mat

    # allocation_list[[i]] <- .pred_alloc <- apply(.mdi_alloc[-(1:burn), ], 2, getmode)
    # compare_df[[dataset_names[i]]] <- .pred_alloc

    compare_tibble$pred_allocation[i + (j - 1) * num_files][[1]] <- .pred_alloc


    # if (i == 1) {
    # row.names(compare_df) <- names(.pred_alloc)
    # }
  }
  # compare_tibble$MDI_allocation[j] <- mdi_allocation
  # compare_tibble$Pred_allocation[j] <- compare_df
}

saveRDS(
  compare_tibble,
  paste0(save_path, "compare_tibble.rds")
)

# === Plot PSM trees ===========================================================
if (do_dendrograms_ie_trees) {
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

# === Phi plots ================================================================

for (i in 1:length(phis)) {
  for (j in 1:ncol(phis[[i]])) {
    curr_phi <- colnames(phis[[i]])[j]
    plot_name <- paste0(save_path, curr_phi, "_plot", plot_type)

    # png(plot_name)
    #
    # phis[[i]][, ..j] %>%
    #   unlist() %>%
    #   plot()
    #
    # dev.off()

    density_plot_name <- paste0(save_path, curr_phi, "_density_plot", plot_type)

    density_title <- paste(curr_phi, ": density plot")

    ggplot(data = phis[[1]], aes_string(x = curr_phi)) +
      geom_density() +
      labs(
        title = density_title
      )

    ggsave(density_plot_name)
  }
}

# ggplot(data = phis[[1]], aes(x = Phi_12)) +
#   geom_density()


# === Check clustering ocnvergence =============================================

# If making heatplots of the clusterings across iterations
if (do_heatplot_clusterings) {

  # Define the type of file to save
  # plot_type <- ".pdf" # ".png" or ".pdf"

  for (j in 1:num_files) {
    # curr_save_path <- paste0(save_path, "seed_", j)
    # Iterate over datasets
    for (i in 1:num_datasets) {
      # Find the dataset name
      dataset <- dataset_names[[i]]

      # Create the save location and file name
      file_name <- paste0(save_path, "similarity_matrix_", dataset, plot_type)

      # Create the title of the plot
      title <- paste("Similarity matrix for", dataset)

      rel_rows <- probes_present_final[, i]

      # re_mat <- .sim_mat
      # re_mat[! rel_rows, ! rel_rows] <- 100
      #
      # ph1 <- pheatmap(re_mat)
      # ord1 <- ph1$tree_col[["order"]]
      # pheatmap(.sim_mat[ord1, ord1], cluster_cols = F, cluster_rows = F)
      #
      # rel_sim_mat <- compare_tibble$similarity_matrix[i + (j - 1) * num_files][[1]] %>%
      #   magrittr::extract(., rel_rows, )

      # Make the heatmap (note that we do not cluster rows).
      # We use a sample of the iterations as otherwise it becomes too heavy
      if (save_plots) {
        ph <- pheatmap(compare_tibble$similarity_matrix[i + (j - 1) * num_files][[1]],
          main = title,
          cluster_rows = T,
          filename = file_name
        )
      } else {
        ph <- pheatmap(compare_tibble$similarity_matrix[i + (j - 1) * num_files][[1]],
          main = title,
          cluster_rows = T
        )
      }
    }
  }
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

# === rand indexing =============================================================
# If instructed to make Rand index plots
if (do_rand_plot) {

  # Define the plot type
  # plot_type <- ".pdf" # ".png" or ".pdf"

  # Create the lists to hold the dataset specific results
  rand <- list()
  rand_plots <- list()


  for (j in 1:num_files) {

    # Loop over the datasets
    for (i in 1:num_datasets) {
      # The current dataset name
      dataset <- dataset_names[[i]]

      # Append this to the generic title to create the specific title
      curr_title <- paste(generic_title, dataset)

      # Similarly for the name of the save file
      curr_save_file <- paste0(generic_save_name, "rand_index_plot_", dataset, plot_type)

      # Create a vector of the adjusted rand index comparing the modal cluster
      # to the clustering at each iteration
      rand[[i + (j - 1) * num_files]] <- apply(
        compare_tibble$mdi_allocation[i + (j - 1) * num_files][[1]],
        # mdi_allocation[[i]],
        1,
        mclust::adjustedRandIndex,
        compare_tibble$mdi_allocation[i][[1]][eff_n_iter, ]
        # compare_df[, i]
      )

      # As we use ggplot2, put this in a dataframe
      plot_data <- data.frame(Rand_index = rand[[i]], Iteration = 1:length(rand[[i]]))

      # Plot the Rand index against iteration number
      rand_plots[[i]] <- ggplot2::ggplot(data = plot_data) +
        ggplot2::geom_point(ggplot2::aes(x = Iteration, y = Rand_index)) +
        ggplot2::labs(
          title = curr_title,
          subtitle = "Comparing modal clustering to clustering at each iteration",
          x = "Index",
          y = "Adjusted Rand Index"
        )

      if (save_plots) {
        # Save the plot
        ggplot2::ggsave(curr_save_file, plot = rand_plots[[i]])
      }
    }
  }
}

if (args$time) {
  print((Sys.time() - stm_i))
}
