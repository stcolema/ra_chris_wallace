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

source("~/Desktop/MDI/mdipp-1.0.1/scripts/analysis.R")
Rcpp::sourceCpp("Analysis/posterior_sim_mat.cpp")
library(magrittr)
library(dplyr)
library(data.table)
library(rlist)
library(pheatmap) # install.packages("pheatmap", dep = T)
library(RColorBrewer)
library(mclust)
library(ggplot2)

# output <- loadDataGauss("~/Desktop/First attempt/output_1.csv")

# === Setup ====================================================================

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
mdi_output_file <- "Analysis/MDI_runs/Full_sets_NAs_dropped_above_0.1/output_full_na_sets.csv"

# Convert this into useable output using functions provided by Sam Mason
mcmc_output <- readMdiMcmcOutput(mdi_output_file)

clust_occ <- getClustersOccupied(mcmc_output)
head(clust_occ)

# By the 8 x 25 = 200 iteration we have reached stationarity in the number of
# clusters occupied
summary(clust_occ[8:800, ])

# Declare empty variable to capture information
n_genes <- NA
mdi_allocation <- list()
allocation_list <- list()

# MDI call specific values
num_datasets <- 9
n_iter <- 2e4
thin <- 25
burn <- 0.1 * (n_iter / thin)

eff_n_iter <- n_iter / thin # - burn

# Create a dataframe to hold the output
col_names <- paste0("D", 1:num_datasets)

# The number of genes is the number of columns in the mcmc_output
# exclusing the columns containing information on the phi and mass parameters
n_genes <- 18523 # SPECIFIC TO CURRENT RUN
if (is.na(n_genes)) {
  n_genes <- mcmc_output %>%
    select(contains("Dataset")) %>%
    ncol() / num_datasets
}

# Create an empty dataframe with column names corresponding to dataset numbers
compare_df <- data.frame(matrix(ncol = num_datasets, nrow = n_genes))
# colnames(compare_df) <- col_names

# The actual names of the datasets are
files_present <- list.files(path = "~/Desktop/MDI/Data/Full_NA_filled_data") %>%
  grep(".csv", ., value = TRUE)

# Remove the output file from this list if you saved it in the same location as
# the data files
if (mdi_output_file %in% files_present) {
  index_to_empty <- files_present == mdi_output_file
  files_present %<>%
    list.remove(index_to_empty)
}
# Remove the file extension
dataset_names <- tools::file_path_sans_ext(files_present)

# Have to have the order from the call
# CD4, CD8, CVD19, CD14, CD15, platelets, ileonic, colonic and rectal biopsies
# transcriptome data for six circulating immune cell types (CD4+ T lymphocytes,
# CD8+ T lymphocytes, CD19+ B lymphocytes, CD14+ monocytes, CD15+ granulocytes,
# platelets) as well as ileal, colonic, and rectal biopsies (IL, TR, RE)
# 323 healthy Europeans
# files_present <- c("CD14.csv", "CD15.csv", "IL.csv", "CD4.csv", "CD8.csv", "CD19.csv", "PLA.csv", "RE.csv", "TR.csv")
dataset_names <- tools::file_path_sans_ext(files_present)

# Set the cell type to the column names (make sure order is as per MDI call)
colnames(compare_df) <- dataset_names

# === MDI output ===============================================================

mdi_pos_sim_mat <- list()
check_median_makes_sense_map <- list()


# x4 <- sample(1:eff_iter, 50, replace=F)

# Capture the allocation information in the named lists and the predicted
# allocation in the dataframe
for (i in 1:num_datasets) {
  dataset_name <- paste0("D", i)

  # Get the allocation and drop the burn in
  mdi_allocation[[i]] <- .mdi_alloc <- getMdiAllocations(mcmc_output, i)

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
  compare_df[[dataset_names[i]]] <- .pred_alloc

  if (i == 1) {
    row.names(compare_df) <- names(.pred_alloc)
  }
}

# === Check clustering ocnvergence =============================================

# If making heatplots of the clusterings across iterations
if (do_heatplot_clusterings) {

  # Define the type of file to save
  # plot_type <- ".pdf" # ".png" or ".pdf"

  # Create the common save path
  save_path <- "Analysis/MDI_runs/Full_sets_NAs_dropped_above_0.1/"

  # Iterate over datasets
  for (i in 1:num_datasets) {
    # Find the dataset name
    dataset <- dataset_names[[i]]

    # Create the save location and file name
    file_name <- paste0(save_path, "iteration_heatplot_before_burn", dataset, plot_type)

    # Create the title of the plot
    title <- paste("Clustering across sample of iterations for", dataset)

    # Make the heatmap (note that we do not cluster rows).
    # We use a sample of the iterations as otherwise it becomes too heavy
    if(save_plots){
    pheatmap(mdi_allocation[[i]][1:80, ],
      main = title,
      cluster_rows = F,
      filename = file_name
    )
    } else {
      pheatmap(mdi_allocation[[i]][c(seq(1, eff_n_iter, 25), eff_n_iter), ],
               main = title,
               cluster_rows = F
      )
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

  # Save the common part of each title
  generic_title <- "MDI: Adjusted Rand index for"
  generic_save_name <- "Analysis/MDI_runs/Full_sets_with_NAs_dropped_if_above_0.1/"

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
    rand[[i]] <- apply(
      mdi_allocation[[i]],
      1,
      mclust::adjustedRandIndex,
      compare_df[, i]
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

    if(save_plots){
      # Save the plot
      ggplot2::ggsave(curr_save_file, plot = rand_plots[[i]])
    }
  }
}

# === Compare individual clusterings to all ====================================

if(do_individual_comparison){
  # Read in the MDI output file
  generic_file_name_start <- "Analysis/MDI_runs/Individual_NAs_dropped_above_0.1/out_transposed_"
  generic_file_end <- "_GE_Corrected4_Covars.csv"
  plot_type <- ".png"
  save_path <- "Analysis/MDI_runs/Individual_NAs_dropped_above_0.1/heatmaps/"
  rand_save_path <- "Analysis/MDI_runs/Individual_NAs_dropped_above_0.1/rand_plots/"
  rand_plots_ind <- list()
  for (i in 1:num_datasets) {
    dataset <- dataset_names[[i]]
    specific_file <- paste0(generic_file_name_start, dataset_names[[i]], generic_file_end)
    specific_mcmc_output <- readMdiMcmcOutput(specific_file)
  
    .mdi_alloc <- getMdiAllocations(specific_mcmc_output, 1) %>%
      .[-(1:burn), ]
  
    # # By checking the imilarity of each row we can decide if the median is an
    # # accurate method to allocate class (if all rows are highly similar label
    # # flipping did not occur)
    .sense_check <- similarity_mat(.mdi_alloc[c(seq(1, eff_n_iter, 25), eff_n_iter), ])
  
    .pred_alloc <- apply(.mdi_alloc, 2, getmode)
  
    names(.pred_alloc)
  
    # Get the clustering from the all MDI
    comp_all_mdi <- compare_df[, i]
    names(comp_all_mdi) <- row.names(compare_df)
  
    order_to_use <- match(names(.pred_alloc), names(comp_all_mdi))
  
  
    compare_clusterings_df <- data.frame(Specific = .pred_alloc, All = comp_all_mdi[order_to_use])
    row.names(compare_clusterings_df) <- names(.pred_alloc)
  
    # Create the save location and file name
    file_name <- paste0(save_path, "mdi_comparison_heatplot", dataset, plot_type)
  
    pheatmap(.sense_check, cluster_rows = F, cluster_cols = F)
    pheatmap(.mdi_alloc[c(seq(1, eff_n_iter, 25), eff_n_iter), ], cluster_rows = F)
    pheatmap_title <- paste0(dataset_names[[i]], ": Comparison of clusterings")
    
    if(save_plots){
      pheatmap(compare_clusterings_df,
        main = dataset_names[[i]],
        filename = file_name
      )
    } else {
      pheatmap(compare_clusterings_df,
               main = dataset_names[[i]]
      )
    }
  
    # Create a vector of the adjusted rand index comparing the modal cluster
    # to the clustering at each iteration
    rand_comp_mode <- apply(
      .mdi_alloc,
      1,
      mclust::adjustedRandIndex,
      .pred_alloc
    )
  
    rand_comp_all <- apply(
      .mdi_alloc,
      1,
      mclust::adjustedRandIndex,
      comp_all_mdi[order_to_use]
    )
  
    # As we use ggplot2, put this in a dataframe
    plot_data_specific <- data.frame(
      Rand_index_self = rand_comp_mode,
      Rand_index_all = rand_comp_all,
      Iteration = 1:length(rand_comp_all)
    )
  
    # Plot the Rand index against iteration number
    rand_plots_ind[[i]] <- ggplot2::ggplot(data = plot_data_specific) +
      ggplot2::geom_point(ggplot2::aes(x = Iteration, y = Rand_index_self)) +
      # ggplot2::geom_point(ggplot2::aes(x = Iteration, y = Rand_index_all)) +
      ggplot2::labs(
        title = dataset,
        subtitle = "Comparing modal clustering to clustering at each iteration",
        x = "Index",
        y = "Adjusted Rand Index"
      )
  
    if(save_plots){
    ggplot2::ggsave(paste0(rand_save_path, dataset, plot_type),
      plot = rand_plots_ind[[i]]
    )
    }
  }
}

# === Some EDA =================================================================

# Possibly need to rearrange based on call order
# compare_df <- compare_df[, c(4,5,1,2,6,7,3,8,9)]

# Check the number and vlaues of clusters present
clusters_present <- apply(compare_df, 2, unique)

# Set the row names
row.names(probes_present_dt) <- probes_present_dt$V1

# Remove the ID column
probes_present_dt_no_V1 <- probes_present_dt[, -1]

# Create a new dataframe with 0's for the probes not present in that dataframe
# and the allocation label otherwise
compare_df_2 <- as.matrix(compare_df) * as.matrix(probes_present_dt_no_V1)

# compare_df_2[compare_df_2 == 0] <- NA

# Add the row and column names
row.names(compare_df_2) <- row.names(compare_df)
colnames(compare_df_2) <- colnames(compare_df)

# Check if this looks correct
head(compare_df_2)

# Move to a new labelling system (as MDI overfits the number of clusters we have
# an unevenly distributed range of values - in our new system there is equal
# distance between each label)
# sort(unique(c(as.matrix(compare_df))))

n_clust <- length(unique(c(compare_df_2)))
old_labels <- sort(unique(c(compare_df_2)))
new_labels <- 1:n_clust
new_labels[1] <- -1
key <- data.frame(old = old_labels, new = new_labels)

# Create a dataframe with the new labels
compare_df_new_labels <- compare_df_2
compare_df_old_labels <- compare_df

# Use [] to ensure that the structure is preserved
compare_df_new_labels[] <- key$new[match(unlist(compare_df_2), key$old)]
compare_df_old_labels[] <- key$new[match(unlist(compare_df), key$old)]

# col_pal <- sort(brewer.pal(15, "Blues"), T)

col_pal_old <- c("#DDDDDD", rev(brewer.pal(n = 11, "RdYlBu"))) # name = "RdYlBu")))
col_pal <- c(rev(brewer.pal(n = 11, "RdYlBu"))) # name = "RdYlBu")))


ph_full <- pheatmap(compare_df_new_labels,
  cluster_rows = T,
  cluster_cols = F,
  color = col_pal_old
)


row_order <- ph_full$tree_row[["order"]]


ph_full_old <- pheatmap(compare_df_old_labels[row_order, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

# Try and find an acceptable colour palette for the heatmap

# col_pal <- rainbow(15, s = 0.6, v = 0.75)
# library(fields) # install.packages("fields")
# col_pal <- tim.colors(15)
# library(gplots) # install.packages("gplots")
# col_pal <- rich.colors(15)

# === Heatmapping ==============================================================

col_pal <- col_pal_old

# Pull out the Gene IDs in the correct order
gene_id <- probe_key %>%
  .[match(.$ProbeID, row.names(compare_df_new_labels))] %>%
  .$Gene

# Move from Probe ID to Gene ID
row.names(compare_df_new_labels) <- gene_id

# Ceck out some specific gene sets
psmd_ind <- grep("PSMD", row.names(compare_df_new_labels))
pheatmap(compare_df_new_labels[psmd_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

PTP4_ind <- grep("PTP4", row.names(compare_df_new_labels))
pheatmap(compare_df_new_labels[PTP4_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

PTPN_ind <- grep("PTPN", row.names(compare_df_new_labels))
pheatmap(compare_df_new_labels[PTPN_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

CFAP_ind <- grep("CFAP", row.names(compare_df_new_labels))
pheatmap(compare_df_new_labels[CFAP_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

NOD_ind <- grep("NOD", row.names(compare_df_new_labels))
pheatmap(compare_df_new_labels[NOD_ind, ],
  luster_rows = F,
  cluster_cols = F,
  color = col_pal
)

ATG_ind <- grep("ATG", row.names(compare_df_new_labels))
pheatmap(compare_df_new_labels[ATG_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

IL_ind <- grep("^IL[1, 2]", row.names(compare_df_new_labels))
pheatmap(compare_df_new_labels[IL_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

MOX_ind <- grep("MOX", row.names(compare_df_new_labels))
pheatmap(compare_df_new_labels[MOX_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

HOX_ind <- grep("HOX", row.names(compare_df_new_labels))
pheatmap(compare_df_new_labels[HOX_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)


# Heatmap of allocaiton
ph_full <- pheatmap(compare_df_new_labels,
  cluster_cols = T,
  color = col_pal
)
row_order <- ph_full$tree_row[["order"]]

pheatmap(compare_df_new_labels[, c(1, 3, 4)], color = col_pal)
pheatmap(compare_df_new_labels[, c(1, 3, 4, 5)], color = col_pal)
pheatmap(compare_df_new_labels[, c(6, 8, 9)], color = col_pal)

df_ph_order <- compare_df_new_labels[row_order, ]

# Inspect this in more manageable section, keeping the same order
# contains(row.names(df_ph_order))
pheatmap(df_ph_order[1:2500, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

pheatmap(df_ph_order[2501:5000, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

# There's very little information here
pheatmap(df_ph_order[5001:10000, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

pheatmap(df_ph_order[10001:15000, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

# Here it quite similar
pheatmap(df_ph_order[15001:18517, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)


pheatmap(compare_df_2, cluster_cols = F)

(n_clusters_present <- lapply(clusters_present, length))

ph_full <- pheatmap(compare_df, cluster_rows = T, cluster_cols = F)
row_order <- ph_full$tree_row[["order"]]
col_order <- ph_full$tree_col[["order"]]

df_ph_order <- compare_df[row_order, ]
pheatmap(df_ph_order[1:2500, ], cluster_rows = F, cluster_cols = F)
pheatmap(df_ph_order[2501:5000, ], cluster_rows = F, cluster_cols = F)
pheatmap(df_ph_order[5001:10000, ], cluster_rows = F, cluster_cols = F)
pheatmap(df_ph_order[10001:15000, ], cluster_rows = F, cluster_cols = F)
pheatmap(df_ph_order[15001:18517, ], cluster_rows = F, cluster_cols = F)

# comparison_sets_1 <- comparison_sets_2 <- dataset_names
#
# compare_df_2 <- compare_df
#
# for (i in 1:num_datasets) {
#   comparison_sets_1 <- comparison_sets_1[-(comparison_sets_1 == dataset_names[i])]
#   comparison_sets_2 <- comparison_sets_2[-(comparison_sets_2 == dataset_names[i])]
#   curr_ds <- dataset_names[[i]]
#   var_1 <- rlang::sym(curr_ds)
#   for (comp_ds in comparison_sets_1) {
#     comparison_sets_2 <- comparison_sets_2[-(comparison_sets_2 == comp_ds)]
#     var_2 <- rlang::sym(comp_ds)
#     for (comp_ds_2 in comparison_sets_2) {
#       var_3 <- rlang::sym(comp_ds_2)
#       new_var <- paste0("Common_", curr_ds, "_", comp_ds, "_", comp_ds_2) %>%
#         rlang::sym()
#       compare_df_2 %<>%
#         mutate(!!new_var := !!var_1 == !!var_2 & !!var_1 == !!var_3)
#     }
#   }
# }
#
# colnames(compare_df_2)
# row.names(compare_df_2) <- row.names(compare_df)
# head(compare_df_2)
# compare_df_3 <- compare_df_2[, !(names(compare_df_2) %in% dataset_names)]
# compare_df_3$All <- apply(compare_df_3, 1, sum)
#
# common_genes <- row.names(compare_df_3)[compare_df_3$All >= 1]
#
#
# # compare_df_2 <- compare_df %>%
# #   mutate(
# #     Common_d1d2 = D1 == D2,
# #     Common_d1d3 = D1 == D3,
# #     Common_d2d3 = D2 == D3,
# #     Common_d1d2d3 = D1 == D2 & D1 == D3
# #   )
#
# com_genes <- compare_df[compare_df_2$Common_d1d2d3 == T, ]
# pheatmap(com_genes)
# row.names(com_genes)
#
# com_all <- sum(compare_df$Common_d1d2d3)
# com_d1d2 <- sum(compare_df$Common_d1d2)
# com_d1d3 <- sum(compare_df$Common_d1d3)
# com_d2d3 <- sum(compare_df$Common_d2d3)
#
# sum_df <- data.frame(
#   All = com_all,
#   D1D2 = com_d1d2,
#   D1D3 = com_d1d3,
#   D2D3 = com_d2d3
# )
#
# # psm_1 <- genPosteriourSimilarityMatrix(mdi_1)
# #
# # # pheatmap(psm_1)
# #
# # consensus_psm_1 <- generateConsensusPSM(mcmc_output)
