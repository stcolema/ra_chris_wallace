
library(tibble)
library(mcclust)
library(magrittr)
library(ggplot2)
library(ggforce)
library(rlist)

true_clusters_2 <- readRDS("~/Desktop/ra_chris_wallace/Yeast/true_clusters_for_yeast_data.rds")

function_dir <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/Analysis_script_functions/"

function_scripts <- c(
  "plot_rand_index.R"
)

for (f in paste0(function_dir, function_scripts)) {
  source(f)
}


comparison_data <- tibble(
  Clustering = list(),
  Adjusted_rand_index = numeric(),
  Dataset = character(),
  Type = character() #,
  # Seed = numeric() # ,
  # Weight = numeric()
)

comparison_data_entry <- comparison_data

types <- c("500", "1000", "long_run")
nice_types <- c("500", "1000", "Long run")
# many_seeds_sub_types <- c("500", "1000")
long_run_seeds <- 1:5

# 1000_weight_1

gen_dir <- "/home/MINTS/sdc56/Desktop/Yeast_3_datasets/MDI_output"

rel_dirs <- list.dirs(gen_dir, recursive = F)
tibble_names <- paste0(rel_dirs, "/compare_tibble.rds")
all_tibbles <- lapply(tibble_names, readRDS)
consensus_types <- c(500, 1000, 5000, 10000)
long_types <- 1:10

datasets <- all_tibbles[[1]]$dataset
types <- c(consensus_types, long_types)

n_types <- length(types)

nice_types <- c(paste0("Consensus ", consensus_types), paste0("Long run: seed ", long_types))

for (j in 1:n_types) {
  curr_type <- types[j]
  curr_nice_type <- nice_types[j]

  curr_tib <- all_tibbles[[j]]
  for (d in datasets) {
    # rel_comp <- comparison_data[which(comparison_data$Dataset == d & comparison_data$Type == "Truth"), ]

    curr_ind <- which(curr_tib$dataset == d)
    
    curr_clustering <- curr_tib$mdi_allocation[[curr_ind]]
    curr_rand_vec <- apply(curr_clustering, 1, unlist_arandi, true_clusters[[d]])
    
    
    comparison_data_entry <- tibble(
      Dataset = d,
      Clustering = list(.curr_clusters),
      Adjusted_rand_index = curr_rand_vec,
      Type = paste0(curr_nice_type) #,
      # Seed = k # ,
      # Weight = curr_weight
    )
    
    comparison_data <- rbind(comparison_data, comparison_data_entry)
  }
}

comparison_data_plot <- comparison_data[comparison_data$Dataset == "MDItestdata1", ]

ggplot(data = comparison_data, aes(x = Type, y = Adjusted_rand_index)) +
  geom_boxplot() +
  # facet_wrap_paginate(~ Dataset, ncol = 1, nrow = 3, page = 1) +
  facet_wrap(~Dataset) +
  theme(axis.text.x = element_text(angle = 30, hjust=1, size=10)) +
  labs(
    title = "Comparison consensus clustering for various iterations to long run (2 million iterations)",
    y = "Adjusted rand index"
    
  )

ggsave("comparison_ARI_yeast_data.png")
