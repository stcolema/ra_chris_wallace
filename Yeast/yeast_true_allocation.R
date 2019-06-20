
library(tibble)
library(mcclust)
library(magrittr)
library(ggplot2)
library(ggforce)

function_dir <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/Analysis_script_functions/"

function_scripts <- c(
  "plot_rand_index.R"
)

for (f in paste0(function_dir, function_scripts)) {
  source(f)
}

my_tib <- readRDS("/home/MINTS/sdc56/Desktop/Yeast_MDI/Yeast_output/Many_seeds_1000_iter/compare_tibble.rds")
weights <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3, 4, 5)

comparison_data <- tibble(
  Clustering = list(),
  Adjusted_rand_index = numeric(),
  Dataset = character(),
  Type = character(),
  Seed = numeric(),
  Weight = numeric()
)

comparison_data_entry <- comparison_data

n_gen <- 10

datasets <- my_tib$dataset
n_datasets <- length(datasets)

clusters <- list()
true_clusters <- list()
for (i in 1:n_datasets) {
  comparison_data_entry_curr <- comparison_data_entry

  curr_dataset <- my_tib$dataset[[i]]

  orig_clusters <- my_tib$pred_allocation[[i]]

  curr_clusters <- unique(unlist(orig_clusters))
  new_cl_1 <- max(curr_clusters) + 1
  new_cl_2 <- max(curr_clusters) + 2

  gen_name <- "New_gene_"

  new_cl_1_names <- paste0(gen_name, 1:n_gen)
  new_cl_2_names <- paste0(gen_name, (1 + n_gen):(2 * n_gen))

  new_clusters <- c(rep(new_cl_1, n_gen), rep(new_cl_2, n_gen)) %>%
    set_names(c(new_cl_1_names, new_cl_2_names))



  true_clusters[[i]] <- .curr_clusters <- c(orig_clusters, new_clusters)

  adj_rand_ind <- unlist_arandi(.curr_clusters, true_clusters[[i]])

  # comparison_data_entry_curr$Dataset[1] <- curr_dataset
  # comparison_data_entry_curr$Clustering[[1]] <- .curr_clusters
  # comparison_data_entry_curr$Adjusted_ran_index[
  #
  comparison_data_entry <- tibble(
    Dataset = curr_dataset,
    Clustering = list(.curr_clusters),
    Adjusted_rand_index = adj_rand_ind,
    Type = "Truth",
    Seed = NA,
    Weight = c(0.0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3, 4, 5)
  )

  comparison_data <- rbind(comparison_data, comparison_data_entry)
}

names(true_clusters) <- my_tib$dataset

weights <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3, 4, 5)

weight_str <- c("0.0", "0.2", "0.4", "0.6", "0.8", "1", "2", "3", "4", "5")

types <- c("500", "1000", "long_run")
# many_seeds_sub_types <- c("500", "1000")
long_run_seeds <- 1:5

# 1000_weight_1

gen_dir <- "/home/MINTS/sdc56/Desktop/Yeast_MDI/Yeast_output/Diffuse/"
tib_name <- "compare_tibble.rds"
gen_many_seeds_name <- "_weight_"
gen_long_run_name <- c("Long_runs/weight_", "_seed_")

n_weights <- length(weights)
n_types <- length(types)


for (j in 1:n_types) {
  curr_type <- types[j]

  # print(j)
  # print(curr_type)

  for (i in 1:n_weights) {
    curr_weight <- weights[i]
    curr_weight_str <- weight_str[i]

    for (d in datasets) {
      rel_comp <- comparison_data[which(comparison_data$Dataset == d & comparison_data$Type == "Truth"), ]

      if (curr_type == "long_run") {
        for (k in long_run_seeds) {
          curr_dir <- paste0(gen_dir, gen_long_run_name[1], curr_weight_str, gen_long_run_name[2], k, "/")
          curr_tib <- readRDS(paste0(curr_dir, tib_name))

          .curr_clusters <- curr_tib$pred_allocation[[which(curr_tib$dataset == d)]]
          adj_rand_ind <- unlist_arandi(.curr_clusters, rel_comp$Clustering[[1]])

          comparison_data_entry <- tibble(
            Dataset = d,
            Clustering = list(.curr_clusters),
            Adjusted_rand_index = adj_rand_ind,
            Type = curr_type, # paste0(curr_type), "_seed_", k),
            Seed = k,
            Weight = curr_weight
          )

          comparison_data <- rbind(comparison_data, comparison_data_entry)
        }
      }

      else {
        curr_dir <- paste0(gen_dir, curr_type, gen_many_seeds_name, curr_weight_str, "/")
        curr_tib <- readRDS(paste0(curr_dir, tib_name))

        .curr_clusters <- curr_tib$pred_allocation[[which(curr_tib$dataset == d)]]

        adj_rand_ind <- unlist_arandi(.curr_clusters, rel_comp$Clustering[[1]])

        comparison_data_entry <- tibble(
          Dataset = d,
          Clustering = list(.curr_clusters),
          Adjusted_rand_index = adj_rand_ind,
          Type = paste0("many_seeds_", curr_type),
          Seed = NA,
          Weight = curr_weight
        )

        comparison_data <- rbind(comparison_data, comparison_data_entry)
      }
    }
  }
}
comparison_data_plot <- comparison_data[comparison_data$Dataset == "MDItestdata1", ]

ggplot(data = comparison_data, aes(x = Type, y = Adjusted_rand_index)) +
  geom_boxplot() +
  # facet_wrap(~Dataset + Weight)
  facet_wrap_paginate(~ Dataset + Weight, ncol = 3, nrow = 4, page = 3)
