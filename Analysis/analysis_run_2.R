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
library(stringr)
# output <- loadDataGauss("~/Desktop/First attempt/output_1.csv")

# Function to plot a heatmap of some subset of genes
gene_subset <- function(data, gene_str, ..., cols_to_use = 1:ncol(data)) {
  # Find the genes of interest within the row names of the data
  genes_of_interest <- grep(gene_str, row.names(data))

  # Plot the heatmap
  ph <- pheatmap(
    compare_mat[genes_of_interest, cols_to_use],
    ...
  )
}


# === Setup ====================================================================

# setwd("~/Desktop/MDI/Data/Full_sets")

compare_df_new_labels <- fread("Data/allocation _data_run_2.csv")

compare_mat <- as.matrix(compare_df_new_labels[, -10])
row.names(compare_mat) <- compare_df_new_labels$V1

# Define the colour palette for the pheatmaps (garish and contrasting)
col_pal <- c("#DDDDDD", rev(brewer.pal(n = 11, "RdYlBu")))

# === Gene Sets ================================================================

# Function to plot a heatmap of some subset of genes
gene_subset <- function(data, gene_str, ..., cols_to_use = 1:ncol(data)) {
  # Find the genes of interest within the row names of the data
  genes_of_interest <- grep(gene_str, row.names(data))

  # Plot the heatmap
  ph <- pheatmap(
    compare_mat[genes_of_interest, cols_to_use],
    ...
  )
}


# Subset the similar datasets
cd_subset <- c(1, 3, 4, 5)
intestine_subset <- c(6, 8, 9)

# Search for the HOX genes
gene_subset(compare_mat, "HOX",
  cluster_rows = T,
  cluster_cols = T,
  color = col_pal,
  main = "HOX genes in CD datasets",
  cols_to_use = cd_subset
)

gene_subset(compare_mat, "HOX",
  cluster_rows = T,
  cluster_cols = T,
  color = col_pal,
  main = "HOX genes in colonic datasets",
  cols_to_use = intestine_subset
)



# Ceck out some specific gene sets
psmd_ind <- grep("PSMD", row.names(compare_mat))
pheatmap(compare_mat[psmd_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

PTP4_ind <- grep("PTP4", row.names(compare_mat))
pheatmap(compare_mat[PTP4_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

PTPN_ind <- grep("PTPN", row.names(compare_mat))
pheatmap(compare_mat[PTPN_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

CFAP_ind <- grep("CFAP", row.names(compare_mat))
pheatmap(compare_mat[CFAP_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

NOD_ind <- grep("NOD", row.names(compare_mat))
pheatmap(compare_mat[NOD_ind, ],
  luster_rows = F,
  cluster_cols = F,
  color = col_pal
)

ATG_ind <- grep("ATG", row.names(compare_mat))
pheatmap(compare_mat[ATG_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

IL_ind <- grep("^IL[1, 2]", row.names(compare_mat))
pheatmap(compare_mat[IL_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

MOX_ind <- grep("MOX", row.names(compare_mat))
pheatmap(compare_mat[MOX_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)

HOX_ind <- grep("HOX", row.names(compare_mat))
pheatmap(compare_mat[HOX_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)


LBP_ind <- grep("LBP", row.names(compare_mat))
pheatmap(compare_mat[LBP_ind, ],
  cluster_rows = F,
  cluster_cols = F,
  color = col_pal
)


# Recognition signals in CD4 according to Wikipedia
TCR_ind <- grep("TCR", row.names(compare_mat))
CD3_ind <- grep("^CD3[A-Z]", row.names(compare_mat))
CD4_ind <- grep("CD4$", row.names(compare_mat))

pheatmap(compare_mat[c(TCR_ind, CD3_ind, CD4_ind), 1:5],
  cluster_rows = T,
  cluster_cols = T,
  color = col_pal
)




# CD31 | PTK7 | CR | IL8
CD_31_ind <- grep("CD31", row.names(compare_mat))
pheatmap(compare_mat[CD_31_ind, ],
  cluster_rows = T,
  cluster_cols = T,
  color = col_pal
)

PTK_ind <- grep("PTK", row.names(compare_mat))
pheatmap(compare_mat[PTK_ind, ],
  cluster_rows = T,
  cluster_cols = T,
  color = col_pal
)

CR_ind <- grep("^CR[1, 2]", row.names(compare_mat))
pheatmap(compare_mat[CR_ind, ],
  cluster_rows = T,
  cluster_cols = T,
  color = col_pal
)

IL_ind <- grep("^IL[5-9]", row.names(compare_mat))
pheatmap(compare_mat[IL_ind, ],
  cluster_rows = T,
  cluster_cols = T,
  color = col_pal
)


# Causal genes

# In the case of IBD, causative genes have been identified for approximately ten
# risk loci on the basis of such “independently” (i.e., not merely reflecting LD
# with other variants) associated coding variants, including:
# NOD2, ATG16L1, IL23R, CARD9, FUT2, and TYK2
# Source: IBD risk loci are enriched in multigenic regulatory modules encompassing putative causative genes
# Yukihide Momozawa et al., 2018. Nature  Communications

NOD_ind <- grep("NOD2", row.names(compare_mat))
ATG_ind <- grep("ATG16L1", row.names(compare_mat))
IL23_ind <- grep("IL23R", row.names(compare_mat))
FUT_ind <- grep("FUT2", row.names(compare_mat))
CARD_ind <- grep("CARD9", row.names(compare_mat))

IBD_genes <- c(NOD_ind, ATG_ind, IL23_ind, FUT_ind, TYK_ind, CARD_ind)
pheatmap(compare_mat[IBD_genes, ],
  cluster_rows = T,
  cluster_cols = T,
  color = col_pal
)

CD_ind <- grep("^CD[1, 4, 8]", row.names(compare_mat))
pheatmap(compare_mat[CD_ind, ],
  cluster_rows = T,
  cluster_cols = T,
  color = col_pal
)


HLA_ind <- grep("HLA", row.names(compare_mat))
pheatmap(compare_mat[HLA_ind, ],
  cluster_rows = T,
  cluster_cols = T,
  color = col_pal
)

IL12_ind <- grep("IL12", row.names(compare_mat))
pheatmap(compare_mat[IL12_ind, ],
  cluster_rows = T,
  cluster_cols = T,
  color = col_pal
)

CCL_ind <- grep("CCL", row.names(compare_mat))
pheatmap(compare_mat[CCL_ind, ],
  cluster_rows = T,
  cluster_cols = T,
  color = col_pal
)

TNF_ind <- grep("TNF", row.names(compare_mat))
pheatmap(compare_mat[TNF_ind, c(1, 3, 4, 5)],
  cluster_rows = T,
  cluster_cols = T,
  color = col_pal
)




gene_subset(compare_mat, "TNF",
  cluster_rows = T,
  cluster_cols = T,
  color = col_pal,
  main = "TNF genes in CD datasets",
  cols_to_use = cd_subset
)

gene_subset(compare_mat, "TNF",
  cluster_rows = T,
  cluster_cols = T,
  color = col_pal,
  main = "TNF genes in colonic datasets",
  cols_to_use = intestine_subset
)



gene_subset(compare_mat, "CCL", cd_subset)
gene_subset(compare_mat, "CCL", intestine_subset)



# === Full dataset visualisation ===============================================

# Heatmap of allocaiton
ph_full <- pheatmap(compare_mat,
  cluster_cols = F,
  color = col_pal
)
row_order <- ph_full$tree_row[["order"]]

df_ph_order <- compare_mat[row_order, ]

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


# Check the number and vlaues of clusters present
clusters_present <- apply(compare_mat, 2, unique)
n_clusters_present <- lapply(clusters_present, length)

# === Similar datasets =========================================================
pheatmap::pheatmap(compare_mat[, c(1, 3, 4, 5)],
  cluster_cols = T,
  color = col_pal,
  main = "Heatmap for CD datasets (all genes)"
)

pheatmap::pheatmap(compare_mat[, c(6, 8, 9)],
  cluster_cols = T,
  color = col_pal,
  main = "Heatmap for intestinal datasets (all genes)"
)

# Less similar

pheatmap::pheatmap(compare_mat[, c(1, 3, 4, 5, 6, 8, 9)],
  cluster_cols = T,
  color = col_pal,
  main = "Heatmap for intestinal & CD datasets (all genes)"
)


pheatmap::pheatmap(compare_mat[, c(7, 2, 1, 3, 4, 5)],
  cluster_cols = T,
  color = col_pal,
  main = "Heatmap for CD datasets + dissimilar datasets"
)

pheatmap::pheatmap(compare_mat[, c(2, 7, 6, 8, 9)],
  cluster_cols = T,
  color = col_pal,
  main = "Heatmap for intestinal datasets + dissimilar datasets"
)

# Pairwise
colnames(compare_mat)
pheatmap::pheatmap(compare_mat[, c(1, 4)],
  cluster_cols = T,
  color = col_pal,
  main = "Heatmap for CD14 and CD4 (all genes)"
)

pheatmap::pheatmap(compare_mat[, c(3, 4)],
  cluster_cols = T,
  color = col_pal,
  main = "Heatmap for CD19 and CD4 (all genes)"
)
