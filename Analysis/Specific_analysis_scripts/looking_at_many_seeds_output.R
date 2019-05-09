#!/usr/env/bin Rscript

# /home/MINTS/sdc56/Desktop/MDI/mdipp-1.0.1/mdipp
# N CD14.csv
# N CD15.csv
# N CD19.csv
# N CD4.csv
# N CD8.csv
# N IL.csv
# N PLA.csv
# RE.csv
# TR.csv
# -n 10000
# -t 25
# -s 1
# > output_full_sets.csv


# Function to find the mode of a vector
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



# === Libraries ================================================================

source("~/Desktop/MDI/mdipp-1.0.1/scripts/analysis.R") # install.packages("mcclust", dep = T)
Rcpp::sourceCpp("Analysis/posterior_sim_mat.cpp") # install.packages("Rcpp", dep = T)
library(magrittr)
library(tibble) # for dataframe of lists
library(dplyr)
library(data.table)
# library(rlist) # install.packages("rlist", dep = T)
library(pheatmap) # install.packages("pheatmap", dep = T)
library(RColorBrewer)
library(mclust) # install.packages("mclust", dep = T)
library(ggplot2)

# output <- loadDataGauss("~/Desktop/First attempt/output_1.csv")

# === Setup ====================================================================

# Have we already done the initial mdi output analysis and prep
data_written <- F

# Instructions to check, save and do several different things
save_plots <- T
plot_type <- ".pdf"
do_heatplot_clusterings <- T
do_rand_plot <- T
do_individual_comparison <- T

# setwd("~/Desktop/MDI/Data/Full_sets")

# Read in the probes present - this is a matrix of bools with the first column
# as the probe IDs and the remaining columns corrresponding to the cell types
# with TRUE indicating the probe is present in this cell type (i.e. not added
# manually with an imputed value) and FALSE indicates we added it in.
probes_present_dt <- fread("Analysis/probes_present_per_dataset.csv")

# Read in the file relating the probe IDs to  the related gene
probe_key <- fread("Analysis/probe_key.csv")

# Read in the MDI output file
file_path <- "/home/MINTS/sdc56/Desktop/Many_seeds_output/"
# file_path <- "Analysis/MDI_runs/vsn_many_seeds/"

# Create the common save path
save_path <- file_path

# For Rand index plots
# Save the common part of each title
generic_title <- "MDI: Adjusted Rand index for"
generic_save_name <- file_path


mdi_output_files <- list.files(path = file_path, full.names = T, include.dirs = F) %>%
  grep("csv", ., value = TRUE)

# mdi_output_files_2 <- "/home/MINTS/sdc56/Desktop/subset_data/Small/Filled_data/output/out_seed_3.csv"
# mdi_output_files <- "/home/MINTS/sdc56/Desktop/matlab_output/CD14_sma_mat_CD19_sma_mat_CD4_sma_mat_CD8_sma_mat_IL_sma_mat_RE_sma_mat_TR_sma_mat_4_mcmcSamples.csv"
# mdi_output_files <- "/home/MINTS/sdc56/Desktop/CD14_small_text_CD4_small_text_2_mcmcSamples.csv"

num_files <- length(mdi_output_files)
file_names <- basename(tools::file_path_sans_ext(mdi_output_files))

# Declare empty variable to capture information
n_genes <- NA
mdi_allocation <- list()
allocation_list <- list()

# MDI call specific values
n_iter <- 1163
thin <- 1
burn <- 0

eff_n_iter <- n_iter / thin # - burn

# For plotting phis
phis <- list()
count <- 0
start_index <- 100 # Consider letting this equal "burn"

# Create a dataframe to hold the output
col_names <- paste0("D", 1:num_datasets)

# The number of genes is the number of columns in the mcmc_output
# exclusing the columns containing information on the phi and mass parameters
n_genes <- 18524 # SPECIFIC TO CURRENT RUN


# The actual names of the datasets are
# files_present <- list.files(path = "~/Desktop/MDI/Data/VSN_NA_data") %>%
# grep(".csv", ., value = TRUE)
files_present <- c(
  "CD14.csv",
  "CD19.csv",
  "CD4.csv",
  "CD8.csv",
  "IL.csv",
  "RE.csv",
  "TR.csv"
)
num_datasets <- length(files_present)

# Remove the file extension
dataset_names <- tools::file_path_sans_ext(files_present)

# mdi_output_file <- "Analysis/MDI_runs/vsn_many_seeds/out_seed_1.csv"
# If already found allocations and written to file, re-load
if (data_written) {
  compare_tibble <- readRDS("~/Desktop/ra_chris_wallace/Analysis/MDI_runs/vsn_many_seeds/compare_tibble.rds")
}


if (!data_written) {
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
    mdi_allocation = rep(vector("list", num_files), num_datasets)
  )

  # === MDI output ===============================================================

  # x4 <- sample(1:eff_iter, 50, replace=F)

  for (j in 1:num_files) {
    # compare_tibble$MDI_allocation[j] <- list()
    # compare_tibble$Pred_allocation[j] <- list()

    mdi_pos_sim_mat <- list()
    check_median_makes_sense_map <- list()

    # Capture the allocation information in the named lists and the predicted
    # allocation in the dataframe


    for (i in 1:num_datasets) {
      dataset_name <- paste0("Dataset", i)

      # Get the allocation and drop the burn in
      # mdi_allocation[[i]] <- .mdi_alloc <- getMdiAllocations(mcmc_out_lst[[j]], i)

      # mdi_alloc2 <- getMdiAllocations(comp_mcmc_out, i)

      mdi_allocation[[i]] <- .mdi_alloc <- mcmc_out_lst[[j]] %>%
        dplyr::select(contains(dataset_name))

      # mdi_pos_sim_mat[[i]][[1]] <- .sim_mat <- similarity_mat(t(.mdi_alloc[,1:1000]))
      mdi_pos_sim_mat[[i]] <- .sim_mat <- similarity_mat(as.matrix(.mdi_alloc))
      pheatmap(.sim_mat)

      if (i == 1) {
        probe_names <- colnames(.mdi_alloc) %>%
          stringr::str_remove_all(paste0(dataset_name, "_")) %>%
          stringr::str_remove_all("X")
      }

      colnames(.mdi_alloc) <- probe_names

      compare_tibble$mdi_allocation[i + (j - 1) * num_files][[1]] <- .mdi_alloc


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

      allocation_list[[i]] <- .pred_alloc <- apply(.mdi_alloc[-(1:burn), ], 2, getmode)
      # compare_df[[dataset_names[i]]] <- .pred_alloc

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
  #   "~/Desktop/ra_chris_wallace/Analysis/MDI_runs/vsn_many_seeds/compare_tibble.rds"
  # )
}
# === Check clustering ocnvergence =============================================

# If making heatplots of the clusterings across iterations
if (do_heatplot_clusterings) {

  # Define the type of file to save
  # plot_type <- ".pdf" # ".png" or ".pdf"

  for (j in 1:num_files) {
    curr_save_path <- paste0(save_path, "seed_", j)
    # Iterate over datasets
    for (i in 1:num_datasets) {
      # Find the dataset name
      dataset <- dataset_names[[i]]

      # Create the save location and file name
      file_name <- paste0(curr_save_path, "iteration_heatplot_before_burn", dataset, plot_type)

      # Create the title of the plot
      title <- paste("Clustering across sample of iterations for", dataset)

      # Make the heatmap (note that we do not cluster rows).
      # We use a sample of the iterations as otherwise it becomes too heavy
      if (save_plots) {
        ph <- pheatmap(compare_tibble$mdi_allocation[i + (j - 1) * num_files][[1]],
          main = title,
          cluster_rows = F,
          filename = file_name
        )
      } else {
        ph <- pheatmap(compare_tibble$mdi_allocation[i + (j - 1) * num_files][[1]],
          main = title,
          cluster_rows = F
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
        compare_tibble$pred_allocation[i + (j - 1) * num_files][[1]]
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

# # === Compare individual clusterings to all ====================================
# 
# if (do_individual_comparison) {
#   # Read in the MDI output file
#   generic_file_name_start <- "Analysis/MDI_runs/Individual_NAs_dropped_above_0.1/out_transposed_"
#   generic_file_end <- "_GE_Corrected4_Covars.csv"
#   plot_type <- ".png"
#   save_path <- "Analysis/MDI_runs/Individual_NAs_dropped_above_0.1/heatmaps/"
#   rand_save_path <- "Analysis/MDI_runs/Individual_NAs_dropped_above_0.1/rand_plots/"
#   rand_plots_ind <- list()
#   for (i in 1:num_datasets) {
#     dataset <- dataset_names[[i]]
#     specific_file <- paste0(generic_file_name_start, dataset_names[[i]], generic_file_end)
#     specific_mcmc_output <- readMdiMcmcOutput(specific_file)
# 
#     .mdi_alloc <- getMdiAllocations(specific_mcmc_output, 1) %>%
#       .[-(1:burn), ]
# 
#     # # By checking the imilarity of each row we can decide if the median is an
#     # # accurate method to allocate class (if all rows are highly similar label
#     # # flipping did not occur)
#     .sense_check <- similarity_mat(.mdi_alloc[c(seq(1, eff_n_iter, 25), eff_n_iter), ])
# 
#     .pred_alloc <- apply(.mdi_alloc, 2, getmode)
# 
#     names(.pred_alloc)
# 
#     # Get the clustering from the all MDI
#     comp_all_mdi <- compare_df[, i]
#     names(comp_all_mdi) <- row.names(compare_df)
# 
#     order_to_use <- match(names(.pred_alloc), names(comp_all_mdi))
# 
# 
#     compare_clusterings_df <- data.frame(Specific = .pred_alloc, All = comp_all_mdi[order_to_use])
#     row.names(compare_clusterings_df) <- names(.pred_alloc)
# 
#     # Create the save location and file name
#     file_name <- paste0(save_path, "mdi_comparison_heatplot", dataset, plot_type)
# 
#     pheatmap(.sense_check, cluster_rows = F, cluster_cols = F)
#     pheatmap(.mdi_alloc[c(seq(1, eff_n_iter, 25), eff_n_iter), ], cluster_rows = F)
#     pheatmap_title <- paste0(dataset_names[[i]], ": Comparison of clusterings")
# 
#     if (save_plots) {
#       pheatmap(compare_clusterings_df,
#         main = dataset_names[[i]],
#         filename = file_name
#       )
#     } else {
#       pheatmap(compare_clusterings_df,
#         main = dataset_names[[i]]
#       )
#     }
# 
#     # Create a vector of the adjusted rand index comparing the modal cluster
#     # to the clustering at each iteration
#     rand_comp_mode <- apply(
#       .mdi_alloc,
#       1,
#       mclust::adjustedRandIndex,
#       .pred_alloc
#     )
# 
#     rand_comp_all <- apply(
#       .mdi_alloc,
#       1,
#       mclust::adjustedRandIndex,
#       comp_all_mdi[order_to_use]
#     )
# 
#     # As we use ggplot2, put this in a dataframe
#     plot_data_specific <- data.frame(
#       Rand_index_self = rand_comp_mode,
#       Rand_index_all = rand_comp_all,
#       Iteration = 1:length(rand_comp_all)
#     )
# 
#     # Plot the Rand index against iteration number
#     rand_plots_ind[[i]] <- ggplot2::ggplot(data = plot_data_specific) +
#       ggplot2::geom_point(ggplot2::aes(x = Iteration, y = Rand_index_self)) +
#       # ggplot2::geom_point(ggplot2::aes(x = Iteration, y = Rand_index_all)) +
#       ggplot2::labs(
#         title = dataset,
#         subtitle = "Comparing modal clustering to clustering at each iteration",
#         x = "Index",
#         y = "Adjusted Rand Index"
#       )
# 
#     if (save_plots) {
#       ggplot2::ggsave(paste0(rand_save_path, dataset, plot_type),
#         plot = rand_plots_ind[[i]]
#       )
#     }
#   }
# }
# 
# # === Some EDA =================================================================
# 
# # Possibly need to rearrange based on call order
# # compare_df <- compare_df[, c(4,5,1,2,6,7,3,8,9)]
# 
# # Check the number and vlaues of clusters present
# clusters_present <- apply(compare_df, 2, unique)
# 
# # Set the row names
# row.names(probes_present_dt) <- probes_present_dt$V1
# 
# # Remove the ID column
# probes_present_dt_no_V1 <- probes_present_dt[, -1]
# 
# # Create a new dataframe with 0's for the probes not present in that dataframe
# # and the allocation label otherwise
# compare_df_2 <- as.matrix(compare_df) * as.matrix(probes_present_dt_no_V1)
# 
# # compare_df_2[compare_df_2 == 0] <- NA
# 
# # Add the row and column names
# row.names(compare_df_2) <- row.names(compare_df)
# colnames(compare_df_2) <- colnames(compare_df)
# 
# # Check if this looks correct
# head(compare_df_2)
# 
# # Move to a new labelling system (as MDI overfits the number of clusters we have
# # an unevenly distributed range of values - in our new system there is equal
# # distance between each label)
# # sort(unique(c(as.matrix(compare_df))))
# 
# n_clust <- length(unique(c(compare_df_2)))
# old_labels <- sort(unique(c(compare_df_2)))
# new_labels <- 1:n_clust
# new_labels[1] <- -1
# key <- data.frame(old = old_labels, new = new_labels)
# 
# # Create a dataframe with the new labels
# compare_df_new_labels <- compare_df_2
# compare_df_old_labels <- compare_df
# 
# # Use [] to ensure that the structure is preserved
# compare_df_new_labels[] <- key$new[match(unlist(compare_df_2), key$old)]
# compare_df_old_labels[] <- key$new[match(unlist(compare_df), key$old)]
# 
# # col_pal <- sort(brewer.pal(15, "Blues"), T)
# 
# col_pal_old <- c("#DDDDDD", rev(brewer.pal(n = 11, "RdYlBu"))) # name = "RdYlBu")))
# col_pal <- c(rev(brewer.pal(n = 11, "RdYlBu"))) # name = "RdYlBu")))
# 
# 
# ph_full <- pheatmap(compare_df_new_labels,
#   cluster_rows = T,
#   cluster_cols = F,
#   color = col_pal_old
# )
# 
# 
# row_order <- ph_full$tree_row[["order"]]
# 
# 
# ph_full_old <- pheatmap(compare_df_old_labels[row_order, ],
#   cluster_rows = F,
#   cluster_cols = F,
#   color = col_pal
# )
# 
# # Try and find an acceptable colour palette for the heatmap
# 
# # col_pal <- rainbow(15, s = 0.6, v = 0.75)
# # library(fields) # install.packages("fields")
# # col_pal <- tim.colors(15)
# # library(gplots) # install.packages("gplots")
# # col_pal <- rich.colors(15)
# 
# # === Compare seeds ============================================================
# 
# row_order <- vector("list", num_datasets)
# cluster_ph_int <- vector("list", num_datasets)
# cluster_ph <- rep(cluster_ph_int, num_files)
# col_pal <- col_pal_old
# 
# # Create the common save path
# save_path <- "Analysis/MDI_runs/vsn_many_seeds/"
# 
# # lie <- matrix(unlist(compare_tibble$Pred_allocation), nrow = 18524, ncol= 7)
# #
# # pheatmap(lie,
# #          color = col_pal,
# #          main = title,
# #          filename = file_name)
# 
# for (i in 1:num_files) {
#   for (j in 1:num_datasets) {
#     curr_save_path <- paste0(save_path, "seed_", j)
# 
#     # Find the dataset name
#     dataset <- dataset_names[[i]]
# 
#     # Create the save location and file name
#     file_name <- paste0(curr_save_path, "predicted_clustering_", dataset, plot_type)
# 
#     # Create the title of the plot
#     title <- paste0(dataset, ": predicted clustering")
# 
#     # If saving plots, save them
#     if (save_plots) {
# 
#       # Take row order from first seed to compare across seeds
#       if (j == 1) {
#         cluster_ph[[i]][[j]] <- pheatmap(compare_tibble$Pred_allocation[[j]],
#           color = col_pal,
#           main = title,
#           filename = file_name
#         )
# 
#         row_order[[i]] <- cluster_ph[[i]][[j]]$tree_row$order
#       } else {
#         cluster_ph[[i]][[j]] <- pheatmap(compare_tibble$Pred_allocation[[j]][row_order[[i]], ],
#           color = col_pal,
#           cluster_rows = F,
#           main = title,
#           filename = file_name
#         )
#       }
#     } else {
#       if (j == 1) {
#         cluster_ph[[i]][[j]] <- pheatmap(compare_tibble$Pred_allocation,
#           color = col_pal
#         )
# 
#         row_order[[i]] <- cluster_ph[[i]][[j]]$tree_row$order
#       } else {
#         cluster_ph[[i]][[j]] <- pheatmap(compare_tibble$Pred_allocation[row_order[[i]], ],
#           color = col_pal,
#           cluster_rows = F
#         )
#       }
#     }
#   }
# }
# 
# # === Heatmapping ==============================================================
# 
# col_pal <- col_pal_old
# 
# # Pull out the Gene IDs in the correct order
# gene_id <- probe_key %>%
#   .[match(.$ProbeID, row.names(compare_df_new_labels))] %>%
#   .$Gene
# 
# # Move from Probe ID to Gene ID
# row.names(compare_df_new_labels) <- gene_id
# 
# # Ceck out some specific gene sets
# psmd_ind <- grep("PSMD", row.names(compare_df_new_labels))
# pheatmap(compare_df_new_labels[psmd_ind, ],
#   cluster_rows = F,
#   cluster_cols = F,
#   color = col_pal
# )
# 
# PTP4_ind <- grep("PTP4", row.names(compare_df_new_labels))
# pheatmap(compare_df_new_labels[PTP4_ind, ],
#   cluster_rows = F,
#   cluster_cols = F,
#   color = col_pal
# )
# 
# PTPN_ind <- grep("PTPN", row.names(compare_df_new_labels))
# pheatmap(compare_df_new_labels[PTPN_ind, ],
#   cluster_rows = F,
#   cluster_cols = F,
#   color = col_pal
# )
# 
# CFAP_ind <- grep("CFAP", row.names(compare_df_new_labels))
# pheatmap(compare_df_new_labels[CFAP_ind, ],
#   cluster_rows = F,
#   cluster_cols = F,
#   color = col_pal
# )
# 
# NOD_ind <- grep("NOD", row.names(compare_df_new_labels))
# pheatmap(compare_df_new_labels[NOD_ind, ],
#   luster_rows = F,
#   cluster_cols = F,
#   color = col_pal
# )
# 
# ATG_ind <- grep("ATG", row.names(compare_df_new_labels))
# pheatmap(compare_df_new_labels[ATG_ind, ],
#   cluster_rows = F,
#   cluster_cols = F,
#   color = col_pal
# )
# 
# IL_ind <- grep("^IL[1, 2]", row.names(compare_df_new_labels))
# pheatmap(compare_df_new_labels[IL_ind, ],
#   cluster_rows = F,
#   cluster_cols = F,
#   color = col_pal
# )
# 
# MOX_ind <- grep("MOX", row.names(compare_df_new_labels))
# pheatmap(compare_df_new_labels[MOX_ind, ],
#   cluster_rows = F,
#   cluster_cols = F,
#   color = col_pal
# )
# 
# HOX_ind <- grep("HOX", row.names(compare_df_new_labels))
# pheatmap(compare_df_new_labels[HOX_ind, ],
#   cluster_rows = F,
#   cluster_cols = F,
#   color = col_pal
# )
# 
# 
# # Heatmap of allocaiton
# ph_full <- pheatmap(compare_df_new_labels,
#   cluster_cols = T,
#   color = col_pal
# )
# row_order <- ph_full$tree_row[["order"]]
# 
# pheatmap(compare_df_new_labels[, c(1, 3, 4)], color = col_pal)
# pheatmap(compare_df_new_labels[, c(1, 3, 4, 5)], color = col_pal)
# pheatmap(compare_df_new_labels[, c(6, 8, 9)], color = col_pal)
# 
# df_ph_order <- compare_df_new_labels[row_order, ]
# 
# # Inspect this in more manageable section, keeping the same order
# # contains(row.names(df_ph_order))
# pheatmap(df_ph_order[1:2500, ],
#   cluster_rows = F,
#   cluster_cols = F,
#   color = col_pal
# )
# 
# pheatmap(df_ph_order[2501:5000, ],
#   cluster_rows = F,
#   cluster_cols = F,
#   color = col_pal
# )
# 
# # There's very little information here
# pheatmap(df_ph_order[5001:10000, ],
#   cluster_rows = F,
#   cluster_cols = F,
#   color = col_pal
# )
# 
# pheatmap(df_ph_order[10001:15000, ],
#   cluster_rows = F,
#   cluster_cols = F,
#   color = col_pal
# )
# 
# # Here it quite similar
# pheatmap(df_ph_order[15001:18517, ],
#   cluster_rows = F,
#   cluster_cols = F,
#   color = col_pal
# )
# 
# 
# pheatmap(compare_df_2, cluster_cols = F)
# 
# (n_clusters_present <- lapply(clusters_present, length))
# 
# ph_full <- pheatmap(compare_df, cluster_rows = T, cluster_cols = F)
# row_order <- ph_full$tree_row[["order"]]
# col_order <- ph_full$tree_col[["order"]]
# 
# df_ph_order <- compare_df[row_order, ]
# pheatmap(df_ph_order[1:2500, ], cluster_rows = F, cluster_cols = F)
# pheatmap(df_ph_order[2501:5000, ], cluster_rows = F, cluster_cols = F)
# pheatmap(df_ph_order[5001:10000, ], cluster_rows = F, cluster_cols = F)
# pheatmap(df_ph_order[10001:15000, ], cluster_rows = F, cluster_cols = F)
# pheatmap(df_ph_order[15001:18517, ], cluster_rows = F, cluster_cols = F)
# 
# # comparison_sets_1 <- comparison_sets_2 <- dataset_names
# #
# # compare_df_2 <- compare_df
# #
# # for (i in 1:num_datasets) {
# #   comparison_sets_1 <- comparison_sets_1[-(comparison_sets_1 == dataset_names[i])]
# #   comparison_sets_2 <- comparison_sets_2[-(comparison_sets_2 == dataset_names[i])]
# #   curr_ds <- dataset_names[[i]]
# #   var_1 <- rlang::sym(curr_ds)
# #   for (comp_ds in comparison_sets_1) {
# #     comparison_sets_2 <- comparison_sets_2[-(comparison_sets_2 == comp_ds)]
# #     var_2 <- rlang::sym(comp_ds)
# #     for (comp_ds_2 in comparison_sets_2) {
# #       var_3 <- rlang::sym(comp_ds_2)
# #       new_var <- paste0("Common_", curr_ds, "_", comp_ds, "_", comp_ds_2) %>%
# #         rlang::sym()
# #       compare_df_2 %<>%
# #         mutate(!!new_var := !!var_1 == !!var_2 & !!var_1 == !!var_3)
# #     }
# #   }
# # }
# #
# # colnames(compare_df_2)
# # row.names(compare_df_2) <- row.names(compare_df)
# # head(compare_df_2)
# # compare_df_3 <- compare_df_2[, !(names(compare_df_2) %in% dataset_names)]
# # compare_df_3$All <- apply(compare_df_3, 1, sum)
# #
# # common_genes <- row.names(compare_df_3)[compare_df_3$All >= 1]
# #
# #
# # # compare_df_2 <- compare_df %>%
# # #   mutate(
# # #     Common_d1d2 = D1 == D2,
# # #     Common_d1d3 = D1 == D3,
# # #     Common_d2d3 = D2 == D3,
# # #     Common_d1d2d3 = D1 == D2 & D1 == D3
# # #   )
# #
# # com_genes <- compare_df[compare_df_2$Common_d1d2d3 == T, ]
# # pheatmap(com_genes)
# # row.names(com_genes)
# #
# # com_all <- sum(compare_df$Common_d1d2d3)
# # com_d1d2 <- sum(compare_df$Common_d1d2)
# # com_d1d3 <- sum(compare_df$Common_d1d3)
# # com_d2d3 <- sum(compare_df$Common_d2d3)
# #
# # sum_df <- data.frame(
# #   All = com_all,
# #   D1D2 = com_d1d2,
# #   D1D3 = com_d1d3,
# #   D2D3 = com_d2d3
# # )
# #
# # # psm_1 <- genPosteriourSimilarityMatrix(mdi_1)
# # #
# # # # pheatmap(psm_1)
# # #
# # # consensus_psm_1 <- generateConsensusPSM(mcmc_output)
