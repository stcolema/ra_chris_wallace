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

setwd("~/Desktop/MDI/Data/Full_sets")

compare_df_new_labels <- fread("allocation _data_run_2.csv")

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

ph_full <- pheatmap(compare_df, cluster_rows = T, cluster_cols = F)
row_order <- ph_full$tree_row[["order"]]
col_order <- ph_full$tree_col[["order"]]

df_ph_order <- compare_df[row_order, ]
pheatmap(df_ph_order[1:2500, ], cluster_rows = F, cluster_cols = F)
pheatmap(df_ph_order[2501:5000, ], cluster_rows = F, cluster_cols = F)
pheatmap(df_ph_order[5001:10000, ], cluster_rows = F, cluster_cols = F)
pheatmap(df_ph_order[10001:15000, ], cluster_rows = F, cluster_cols = F)
pheatmap(df_ph_order[15001:18517, ], cluster_rows = F, cluster_cols = F)
