#!/usr/bin/env Rscript

# Function to plot heatmaps for fused genes from permutations of datasets

fused_gene_heatmaps <- function(expression_data_lst,
                                fused_probes_lst,
                                gene_id,
                                datasets,
                                file_path,
                                n_datasets,
                                plot_type = ".png") {

  # Set up the directory to save too
  dir_name <- paste0(file_path, "Fusion_expression_data/")
  dir.create(dir_name, showWarnings = FALSE)

  # The generic heatmap file name
  generic_ph_title <- paste0(dir_name, "heatmap_")

  # Do the various combinations of datasets
  for (i in 1:(n_datasets - 1)) {
    d_i <- datasets[[i]]
    expression_data_i <- expression_data_lst[[i]]

    for (j in (i + 1):n_datasets) {
      d_j <- datasets[[j]]
      expression_data_j <- expression_data_lst[[j]]

      # Extract the indices for the "fused" and "unfused" genes
      fused_ind <- fused_probes_lst[[i]][[j]] # fused_non_zero_probes[[1]][[2]]
      non_fused_ind <- !fused_ind


      # Find a nice ordering of the columns. This is superfluous and inefficient
      # Also, this seems to be saving as a pdf. Weird.
      # ph1 <- pheatmap(expression_data_i, cluster_rows = F)
      # ph2 <- pheatmap(expression_data_j, cluster_rows = F)
      #
      # col_order_1 <- ph1$tree_col$order
      # col_order_2 <- ph2$tree_col$order
      #
      # new_expression_data <- cbind(
      #   expression_data_i[, col_order_1],
      #   expression_data_j[, col_order_2]
      # ) %>%
      #   magrittr::set_rownames(gene_id)

      # We ignore the above and instead have data unordered
      col_order_1 <- colnames(expression_data_i) # to allow no edits and to uncomment above for seamless integration

      new_expression_data <- cbind(
        expression_data_i,
        expression_data_j
      ) %>%
        magrittr::set_rownames(gene_id)

      # Filename for fused genes heatmap
      fused_ph_file_name <- paste0(
        generic_ph_title,
        "fused_genes_",
        d_i,
        "_",
        d_j,
        plot_type
      )

      # If more than 1 fused gene can attempt to cluster rows
      if (sum(fused_ind) > 1) {
        pheatmap(new_expression_data[fused_ind, ],
          cluster_cols = F,
          gaps_col = length(col_order_1),
          main = paste("Fused probes for", d_i, "and", d_j),
          filename = fused_ph_file_name
        )
      } else {
        # If exactly 1 fused gene can still heatmap but no clustering
        if (sum(fused_ind) == 1) {
          pheatmap(new_expression_data[fused_ind, ],
            cluster_rows = F,
            cluster_cols = F,
            gaps_col = length(col_order_1),
            main = paste("Fused probes for", d_i, "and", d_j),
            filename = fused_ph_file_name
          )
        }
      }

      # File name for the heatmap fo the expression data for the unfused genes
      unfused_ph_file_name <- paste0(
        generic_ph_title,
        "unfused_genes_",
        d_i,
        "_",
        d_j,
        plot_type
      )

      # If have n > 1 can cluster rows. If n >= 1 can heatmap.
      if (sum(non_fused_ind) > 1) {
        pheatmap(new_expression_data[non_fused_ind, ],
          cluster_cols = F,
          gaps_col = length(col_order_1),
          main = paste("Unfused probes for", d_i, "and", d_j),
          filename = unfused_ph_file_name
        )
      } else {
        if (sum(non_fused_ind) == 1) {
          pheatmap(new_expression_data[non_fused_ind, ],
            cluster_cols = F,
            cluster_rows = F,
            gaps_col = length(col_order_1),
            main = paste("Unfused probes for", d_i, "and", d_j),
            filename = unfused_ph_file_name
          )
        }
      }
    }
  }
}
