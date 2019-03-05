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

# === Libraries ================================================================

source("~/Desktop/MDI/mdipp-1.0.1/scripts/analysis.R")
library(magrittr)
library(dplyr)
library(data.table)
library(rlist)
library(pheatmap) # install.packages("pheatmap", dep = T)
library(RColorBrewer)
# output <- loadDataGauss("~/Desktop/First attempt/output_1.csv")

# === Setup ====================================================================

setwd("~/Desktop/MDI/Runs/Run 1")

# Read in the probes present - this is a matrix of bools with the first column
# as the probe IDs and the remaining columns corrresponding to the cell types
# with TRUE indicating the probe is present in this cell type (i.e. not added
# manually with an imputed value) and FALSE indicates we added it in.
probes_present_dt <- fread("~/Desktop/MDI/Data/probes_present_per_dataset.csv")

# Read in the file relating the probe IDs to  the related gene
probe_key <- fread("~/Desktop/MDI/Data/probe_key.csv")

# Read in the MDI output file
mdi_output_file <- "output_all.csv"

# Convert this into useable output using functions provided by Sam Mason
mcmc_output <- readMdiMcmcOutput(mdi_output_file)

# Declare empty variable to capture information
n_genes <- NA
mdi_allocation <- list()
allocation_list <- list()

# MDI call specific values
num_datasets <- 9
n_iter <- 1e4
thin <- 25
burn <- 0.1 * (n_iter / thin)
eff_n_iter <- n_iter / thin - burn

# Create a dataframe to hold the output
col_names <- paste0("D", 1:num_datasets)

# The number of genes is the number of columns in the mcmc_output
# exclusing the columns containing information on the phi and mass parameters
if (is.na(n_genes)) {
  n_genes <- mcmc_output %>%
    select(contains("Dataset")) %>%
    ncol() / num_datasets
}

# Create an empty dataframe with column names corresponding to dataset numbers
compare_df <- data.frame(matrix(ncol = num_datasets, nrow = n_genes))
# colnames(compare_df) <- col_names

# The actual names of the datasets are
files_present <- list.files(path = "~/Desktop/MDI/Data/Prepped_data") %>%
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
files_present <- c("CD14.csv", 
                   "CD15.csv", 
                   "IL.csv",
                   "CD4.csv", 
                   "CD8.csv",
                   "CD19.csv",
                   "PLA.csv",
                   "RE.csv",
                   "TR.csv")

dataset_names <- tools::file_path_sans_ext(files_present)

# Set the cell type to the column names (make sure order is as per MDI call)
colnames(compare_df) <- dataset_names

# === MDI output ===============================================================

# Capture the allocation information in the named lists and the predicted
# allocation in the dataframe
for (i in 1:num_datasets) {
  dataset_name <- paste0("D", i)

  # Get the allocation and drop the burn in
  mdi_allocation[[i]] <- .mdi_alloc <- getMdiAllocations(mcmc_output, i) %>%
    .[-(1:burn), ]

  allocation_list[[i]] <- .pred_alloc <- apply(.mdi_alloc, 2, median)

  compare_df[[dataset_names[i]]] <- .pred_alloc

  if (i == 1) {
    row.names(compare_df) <- names(.pred_alloc)
  }
}

# Possibly need to rearrange based on call order
# compare_df <- compare_df[, c(4,5,1,2,6,7,3,8,9)]

# Check the number and vlaues of clusters present
clusters_present <- apply(compare_df, 2, unique)

n_clust <- length(unique(c(as.matrix(compare_df))))
old_labels <- sort(unique(c(as.matrix(compare_df))))
new_labels <- 1:n_clust
key <- data.frame(old = old_labels, new = new_labels)

# Create a dataframe with the new labels
compare_df_new_labels <- as.matrix(compare_df)
# compare_df_old_labels <- compare_df

# Use [] to ensure that the structure is preserved
compare_df_new_labels[] <- key$new[match(unlist(compare_df), key$old)]
# compare_df_old_labels[] <- key$new[match(unlist(compare_df), key$old)]

# col_pal <- sort(brewer.pal(15, "Blues"), T)

# Write the data
compare_df_new_labels_write <- data.table(compare_df_new_labels)
compare_df_new_labels_write$V1 <- row.names(compare_df_new_labels)
  
fwrite(compare_df_new_labels_write, "allocation_data_run_1.csv")

col_pal <- c("#DDDDDD", rev(brewer.pal(n = 11, "RdYlBu"))) # name = "RdYlBu")))

ph_full <- pheatmap(compare_df_new_labels, 
                    cluster_cols = F,
                    color = col_pal,
                    main = "Heatmap for common probes"
                    )

# === Heatmapping ==============================================================

# Pull out the Gene IDs in the correct order
gene_id <- probe_key %>%
  .[match(.$ProbeID, row.names(compare_df_new_labels))] %>%
  .$Gene %>% 
  na.omit()

# Move from Probe ID to Gene ID
row.names(compare_df_new_labels) <- gene_id
# gene_id <- c(row.names(compare_df_new_labels))

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

IL_ind <- grep("^IL", row.names(compare_df_new_labels))
pheatmap(compare_df_new_labels[IL_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

# Heatmap of allocaiton
ph_full <- pheatmap(compare_df_new_labels,
  cluster_cols = F,
  color = col_pal
)
row_order <- ph_full$tree_row[["order"]]

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

ph_full <- pheatmap(compare_df_new_labels, cluster_rows = T, cluster_cols = F)
row_order <- ph_full$tree_row[["order"]]
col_order <- ph_full$tree_col[["order"]]

df_ph_order <- compare_df_new_labels[row_order, ]
pheatmap(df_ph_order[1:1000, ], cluster_rows = F, cluster_cols = F)
pheatmap(df_ph_order[1001:2000, ], cluster_rows = F, cluster_cols = F)
pheatmap(df_ph_order[2001:3000, ], cluster_rows = F, cluster_cols = F)
pheatmap(df_ph_order[3001:4000, ], cluster_rows = F, cluster_cols = F)
pheatmap(df_ph_order[4001:4964, ], cluster_rows = F, cluster_cols = F)

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
