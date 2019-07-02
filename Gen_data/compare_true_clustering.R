
library(mcclust)
library(magrittr)
library(ggplot2)
library(ggforce)


# For tibbles
library(tibble) # for dataframe of lists

# For data wrangling
library(dplyr)

# Heatmapping
library(pheatmap) # install.packages("pheatmap", dep = T)

# Colour palettes
library(RColorBrewer)


function_dir <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/Analysis_script_functions/"

function_scripts <- c(
  "plot_rand_index.R"
)

for (f in paste0(function_dir, function_scripts)) {
  source(f)
}

gen_dir <- "/home/MINTS/sdc56/Desktop/Gen_data_output"

all_dirs_in_full <- list.dirs(path = gen_dir, recursive = F)

longer_dirs <- all_dirs_in_full %>% 
  extract(grep("Longer" , .))

consensus_dirs <- all_dirs_in_full %>% 
  extract(grep("Many" , .)) 

# gen_long_run_name <- c("Long", "Longer")
# gen_many_seeds_name <- "Many_seeds_"


all_rel_dirs <- c(longer_dirs, consensus_dirs)
n_long <- length(longer_dirs)
long_types <- rep("Long", n_long)
long_types <- paste("Long", 1:n_long)
consensus_types <- paste("Consensus", c(1000, 10000, 500, 5000))
types <- c(long_types, consensus_types)

tib_name <- "compare_tibble.rds"
all_data <- paste0(all_rel_dirs, "/", tib_name)

all_tibs <- lapply(all_data, readRDS)

# === Set up to record =========================================================

# my_tib <- readRDS("/home/MINTS/sdc56/Desktop/Yeast_MDI/Yeast_output/Many_seeds_1000_iter/compare_tibble.rds")
# weights <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3, 4, 5)

true_labels <- c(rep(1, 25),
                 rep(2, 50),
                 rep(3, 75),
                 rep(4, 100),
                 rep(5, 150)
)

comparison_data <- tibble(
  # Clustering = list(),
  Adjusted_rand_index = numeric(),
  Dataset = character(),
  Type = character()
  # Seed = numeric()
)

comparison_data_entry <- comparison_data

# 1000_weight_1

datasets <- all_tibs[[1]]$dataset
n_datasets <- length(datasets)

n_mdi <- length(all_data)
for(i in 1:n_mdi){
  curr_tib <- all_tibs[[i]]
  curr_type <- types[i]
for(d in datasets){
  curr_index <- which(curr_tib$dataset == d)
  
  curr_clustering <- curr_tib$mdi_allocation[[curr_index]]
  
  # curr_clustering <- curr_tib$pred_allocation[[curr_index]] %>% unname()
  # curr_adj_rand <- mcclust::arandi(true_labels, curr_clustering)
  
  curr_rand_vec <- apply(curr_clustering, 1, unlist_arandi, true_labels)
  
  
  comparison_data_entry <- tibble(
    Adjusted_rand_index = curr_rand_vec,
    Dataset = d,
    Type = curr_type
  )
  
  # cat("\n\nDataset: ", d)
  # cat("\nType: ", curr_type)
  # cat("\nRand index: ", curr_adj_rand)
  # cat("\nCurr clustering (head)\n: ", head(curr_clustering))
  # cat("\n\nCurr clustering (str)\n: ", str(curr_clustering))
  
  comparison_data <- rbind(comparison_data, comparison_data_entry)
}
}

# small_subset <- c(11:14)
# small_tibs <- all_tibs[small_subset]
# small_types <- types[small_subset]
# small_comparison <-  tibble(
#   Adjusted_rand_index = numeric(),
#   Dataset = character(),
#   Type = character()
# )
# 
# n_mdi_small <- length(small_types)
# for(i in 1:n_mdi_small){
#   curr_tib <- small_tibs[[i]]
#   curr_type <- small_types[i]
#   for(d in datasets){
#     curr_index <- which(curr_tib$dataset == d)
#     curr_clustering <- curr_tib$mdi_allocation[[curr_index]]
#     curr_rand_vec <- apply(curr_clustering, 1, unlist_arandi, true_labels)
# 
#     # curr_adj_rand <- mcclust::arandi(true_labels, curr_clustering)
#     comparison_data_entry <- tibble(
#       Adjusted_rand_index = curr_rand_vec,
#       Dataset = d,
#       Type = curr_type
#     )
#     
#     small_comparison <- rbind(small_comparison, comparison_data_entry)
#   }
# }

# mdi_alloc <- all_tibs[[1]]$mdi_allocation[[1]][1:10, 1:10]
# 
# 
# rand[[i + (j - 1) * n_files]] <- apply(
#   mdi_allocation_lst[i + (j - 1) * n_files][[1]],
#   # mdi_allocation[[i]],
#   1,
#   # mcclust::arandi,
#   unlist_arandi,
#   # mclust::adjustedRandIndex,
#   mdi_allocation_lst[i][[1]][eff_n_iter - start_index, ]
#   # compare_df[, i]
# )
# 
# 
# consesnsus_tibble <- comparison_data[comparison_data$Type %in% consensus_types, ]

ggplot(data = comparison_data, aes(x = Type, y = Adjusted_rand_index)) +
  geom_violin() +
  facet_wrap(~Dataset) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  labs(
    title = "Comparison consensus clustering for various iterations to long run (2 million iterations)",
    y = "Adjusted rand index"
    
  )
  # facet_wrap_paginate(~ Dataset + Weight)

comparison_data_2 <- comparison_data

comparison_data_2$Type[comparison_data_2$Type %in% long_types] <- "Long"

ggplot(data = comparison_data, aes(x = Type, y = Adjusted_rand_index)) +
  geom_boxplot(notch = T) +
  facet_wrap(~Dataset) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  labs(
    title = "Comparison consensus clustering for various iterations to long run (2 million iterations)",
    y = "Adjusted rand index"
    
  )

ggplot(data = comparison_data, aes(x = Type, y = Adjusted_rand_index)) +
  geom_boxplot() +
  facet_wrap_paginate(~Dataset, nrow = 1, ncol = 1) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  labs(
    title = "Comparison consensus clustering for various iterations to long run (2 million iterations)",
    y = "Adjusted rand index"
    
  )

ggplot(data = comparison_data_2, aes(x = Type, y = Adjusted_rand_index)) +
  geom_violin() +
  facet_wrap(~Dataset) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  labs(
    title = "Comparison consensus clustering for various iterations to long run (2 million iterations)",
    y = "Adjusted rand index"
    
  )


# 
# ggplot(data = small_comparison, aes(x = Type, y = Adjusted_rand_index)) +
#   geom_boxplot() +
#   facet_wrap(~Dataset) +
#   theme(axis.text.x = element_text(angle = 30, hjust=1)) +
#   labs(
#     title = "Comparison consensus clustering for various iterations to long run (2 million iterations)",
#     y = "Adjusted rand index"
#     
#   )
# 
# 
# ggplot(data = consesnsus_tibble, aes(x = Type, y = Adjusted_rand_index)) +
#   geom_point() +
#   facet_wrap(~Dataset) +
#   theme(axis.text.x = element_text(angle = 30, hjust=1)) +
#   labs(
#     title = "Comparison consensus clustering for various iterations to long run (2 million iterations)",
#     y = "Adjusted rand index"
#     
#   )
