


function_dir <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/Analysis_script_functions/"

function_scripts <- c(
  "plot_similarity_matrices.R",
  "plot_comparison_expression_clustering.R",
  "create_psm.R"
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

yeast_dir <- "/home/MINTS/sdc56/Desktop/Yeast_MDI/Yeast_output/Diffuse/"
gen_long_runs_dir <- c("Long_runs/weight_", "_seed_")
gen_many_seeds_dir <- "_weight_"

tib_name <- "compare_tibble.rds"

weights <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3, 4, 5)

weight_str <- c("0.0", "0.2", "0.4", "0.6", "0.8", "1", "2", "3", "4", "5")

types <- c("500", "1000", "long_run")
# many_seeds_sub_types <- c("500", "1000")
long_run_seeds <- 1:5

n_weights <- length(weights)
n_types <- length(types)

datasets <- c("MDItestdata1", "MDItestdata2", "MDItestdata3", "MDItestdata4", "MDItestdata5", "MDItestdata6")

col_pal_sim <- colorRampPalette(c("white", "#146EB4"))(100)
palette_length <- length(col_pal)

sim_breaks <- c(
  seq(1 / palette_length, 1, length.out = palette_length)
)

col_pal_expr <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(100)

expr_breaks <- define_breaks(col_pal_expr)

sim_comparison_tibble <- tibble(
  Dataset = character(),
  Similarity_matrix = list(),
  Correlation_matrix = list(),
  Expression_data = list(),
  Similarity_matrix_ph = list(),
  Correlation_matrix_ph = list(),
  Expression_data_ph = list(),
  Sim_row_order = list(),
  Expr_col_order = list(),
  Weight = numeric(),
  Type = character(),
  Seed = numeric()
)



for (j in 1:n_types) {
  curr_type <- types[j]

  for (i in 1:n_weights) {
    curr_weight <- weights[i]
    curr_weight_str <- weight_str[i]


    if (curr_type == "long_run") {
      for (k in long_run_seeds) {
        curr_dir <- paste0(gen_dir, gen_long_run_name[1], curr_weight_str, gen_long_run_name[2], k, "/")
        curr_tib <- readRDS(paste0(curr_dir, tib_name))

        for (d in datasets) {
          .curr_index <- which(curr_tib$dataset == d)

          .curr_sim <- curr_tib$similarity_matrix[[.curr_index]]

          row_order <- sim_comparison_tibble$Sim_row_order[[which(sim_comparison_tibble$Dataset == d
          & sim_comparison_tibble$Type == "500"
          & sim_comparison_tibble$Weight == curr_weight)]]
          expr_col_order <- sim_comparison_tibble$Expr_col_order[[which(sim_comparison_tibble$Dataset == d
          & sim_comparison_tibble$Type == "500"
          & sim_comparison_tibble$Weight == curr_weight)]]


          .curr_sim <- .curr_sim[row_order, row_order]
          .curr_corr <- curr_tib$correlation_matrix[[.curr_index]][row_order, row_order]
          .curr_expr <- curr_tib$expression_data[[.curr_index]][row_order, expr_col_order]

          ph_sim <- pheatmap(.curr_sim,
            cluster_rows = F,
            cluster_cols = F,
            color = col_pal_sim,
            breaks = sim_breaks
          )

          ph_corr <- pheatmap(.curr_corr,
            cluster_rows = F,
            cluster_cols = F,
            color = col_pal_expr,
            breaks = expr_breaks
          )

          ph_expr <- pheatmap(.curr_expr,
            cluster_rows = F,
            cluster_cols = F,
            color = col_pal_expr,
            breaks = expr_breaks
          )

          # Add to a tibble
          sim_comparison_tibble_entry <- tibble(
            Dataset = d,
            Similarity_matrix = list(.curr_sim),
            Correlation_matrix = list(.curr_corr),
            Expression_data = list(.curr_expr),
            Similarity_matrix_ph = list(ph_sim),
            Correlation_matrix_ph = list(ph_corr),
            Expression_data_ph = list(ph_expr),
            Sim_row_order = list(row_order),
            Expr_col_order = list(expr_col_order),
            Weight = curr_weight,
            Type = curr_type,
            Seed = k
          )

          sim_comparison_tibble <- rbind(sim_comparison_tibble, sim_comparison_tibble_entry)
        }
      }
    }

    else {

      # Find the current directory
      curr_dir <- paste0(gen_dir, curr_type, gen_many_seeds_name, curr_weight_str, "/")

      # Read in the tibble containing the data relevant to the current format
      curr_tib <- readRDS(paste0(curr_dir, tib_name))

      for (d in datasets) {

        # The relevant index within said tibble
        .curr_index <- which(curr_tib$dataset == d)

        # Extract the similarity matrix and expression data
        .curr_sim <- curr_tib$similarity_matrix[[.curr_index]]
        .curr_expr <- curr_tib$expression_data[[.curr_index]]

        if (curr_type == "500") {
          # Find the row order based on the similarity matrix
          sim_hc <- hclust(dist(.curr_sim))
          row_order <- sim_hc$order

          # Find the column order for the expression data
          expr_hc <- hclust(dist(t(toy_expr)))
          expr_col_order <- expr_hc$order
        } else {
          row_order <- sim_comparison_tibble$Sim_row_order[[which(sim_comparison_tibble$Dataset == d
          & sim_comparison_tibble$Type == "500"
          & sim_comparison_tibble$Weight == curr_weight)]]
          expr_col_order <- sim_comparison_tibble$Expr_col_order[[which(sim_comparison_tibble$Dataset == d
          & sim_comparison_tibble$Type == "500"
          & sim_comparison_tibble$Weight == curr_weight)]]
        }

        # Rearragne data
        .curr_sim <- .curr_sim[row_order, row_order]
        .curr_corr <- curr_tib$correlation_matrix[[.curr_index]][row_order, row_order]
        .curr_expr <- .curr_expr[row_order, expr_col_order]

        # Create heatmaps
        ph_sim <- pheatmap(.curr_sim,
          cluster_rows = F,
          cluster_cols = F,
          color = col_pal_sim,
          breaks = sim_breaks
        )

        ph_corr <- pheatmap(.curr_corr,
          cluster_rows = F,
          cluster_cols = F,
          color = col_pal_expr,
          breaks = expr_breaks
        )

        ph_expr <- pheatmap(.curr_expr,
          cluster_rows = F,
          cluster_cols = F,
          color = col_pal_expr,
          breaks = expr_breaks
        )

        # Add to a tibble
        sim_comparison_tibble_entry <- tibble(
          Dataset = d,
          Similarity_matrix = list(.curr_sim),
          Correlation_matrix = list(.curr_corr),
          Expression_data = list(.curr_expr),
          Similarity_matrix_ph = list(ph_sim),
          Correlation_matrix_ph = list(ph_corr),
          Expression_data_ph = list(ph_expr),
          Sim_row_order = list(row_order),
          Expr_col_order = list(expr_col_order),
          Weight = curr_weight,
          Type = curr_type,
          Seed = NA
        )

        # Save to big tibble
        sim_comparison_tibble <- rbind(sim_comparison_tibble, sim_comparison_tibble_entry)
      }
    }
  }
}

curr_dataset <- datasets[1]
curr_weight <- weights[3]

ph_list_1 <- sim_comparison_tibble$Similarity_matrix_ph[which(sim_comparison_tibble$Dataset == curr_dataset
& sim_comparison_tibble$Weight == curr_weight)]

n_ph <- length(ph_list_1)

curr_corr_ph <- sim_comparison_tibble$Correlation_matrix_ph[which(sim_comparison_tibble$Dataset == curr_dataset
& sim_comparison_tibble$Weight == curr_weight)][[1]]

for (i in 2:n_ph) {
  save_name <- paste0(
    "/home/MINTS/sdc56/Desktop/ph_comp_",
    curr_dataset,
    "_weight_",
    curr_weight,
    "_500_many_",
    i,
    ".png"
  )

  main_title <- paste0(
    "Comparison of similarity matrices for ",
    curr_dataset,
    " for weight ",
    curr_weight,
    " (correlation matrix on right)"
  )

  combine_pheatmaps(list(ph_list_1[[1]]$gtable, ph_list_1[[i]]$gtable, curr_corr_ph$gtable),
    save_name = save_name,
    main = main_title
  )
}
