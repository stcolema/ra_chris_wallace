

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

# To load true clustering
library(rlist)


function_dir <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/Analysis_script_functions/"

function_scripts <- c(
  "plot_rand_index.R"
)

for (f in paste0(function_dir, function_scripts)) {
  source(f)
}

# For controlling ordering of boxplots
add_level <- function(x, new_level){
  if(is.factor(x)) return(factor(x, levels=c(levels(x), new_level)))
  return(x)
}

# === Main script ==============================================================

# Set ggplot2 theme
theme_set(theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust=1)))

true_clusters <- readRDS("~/Desktop/ra_chris_wallace/Yeast/true_clusters_for_yeast_data.rds")



comparison_data <- tibble(
  # Clustering = list(),
  Adjusted_rand_index = numeric(),
  Dataset = character(),
  Type = character() #,
  # Seed = numeric() # ,
  # Weight = numeric()
)

comparison_data_entry <- comparison_data

# types <- c("500", "1000", "long_run")
# nice_types <- c("500", "1000", "Long run")
# many_seeds_sub_types <- c("500", "1000")
long_run_seeds <- 1:5

# 1000_weight_1

gen_dir <- "/home/MINTS/sdc56/Desktop/Yeast_3_datasets/MDI_output"

all_dirs_in_full <- list.dirs(path = gen_dir, recursive = F)
# all_dirs_short <- list.dirs(path = gen_dir, recursive = F, full.names = F)

longer_dirs <- all_dirs_in_full %>% 
  extract(grep("Long" , .))

consensus_dirs <- all_dirs_in_full %>% 
  extract(grep("Consensus" , .)) 

tib_types <- list.dirs(path = gen_dir, recursive = F, full.names = F)

# gen_long_run_name <- c("Long", "Longer")
# gen_many_seeds_name <- "Many_seeds_"


all_rel_dirs <- all_dirs_in_full # c(consensus_dirs, longer_dirs)
# n_long <- length(longer_dirs)
# long_types <- rep("Long", n_long)
# long_types <- paste("Long", 1:n_long)
# consensus_types <- paste("Consensus", c(1000, 10000, 500, 5000))
consensus_types <- all_dirs_short %>% 
  extract(grep("Consensus" , .)) 

longer_dirs <- all_dirs_short %>% 
  extract(grep("Long" , .))


# types <- c(consensus_types, long_types)

types <- stringr::str_replace(tib_types, "_", " ")
n_types <- length(types)

# types <- types[match(stringr::str_replace(tib_types, "_", " "), types)]

tib_name <- "compare_tibble.rds"
all_data <- paste0(all_rel_dirs, "/", tib_name)

all_tibs <- lapply(all_data, readRDS)

types <- tib_types

consensus_types <- tib_types %>% extract(grep("Consensus", .))
long_types <- tib_types %>% extract(grep("Long", .))

long_nice_types <- str_replace(long_types, "_", " run: seed ")
nice_types <- c(str_replace(consensus_types, "_", " "), str_replace(long_types, "_", " run: seed " ))

# === Recording data ===========================================================

for (j in 1:n_types) {
  curr_type <- types[j]
  curr_nice_type <- nice_types[j]

  curr_tib <- all_tibs[[j]]
  for (d in datasets) {
    # rel_comp <- comparison_data[which(comparison_data$Dataset == d & comparison_data$Type == "Truth"), ]

    curr_ind <- which(curr_tib$dataset == d)
    
    curr_clustering <- curr_tib$mdi_allocation[[curr_ind]]
    
    if(grepl("Long", curr_type)){
      burn_in <- 10
      n <- nrow(curr_clustering)
      curr_rand_vec <- apply(curr_clustering[burn_in:n, ], 1, unlist_arandi, true_clusters[[d]])
      
    } else {
      curr_rand_vec <- apply(curr_clustering, 1, unlist_arandi, true_clusters[[d]])
    }
    # curr_rand_vec <- apply(curr_clustering, 1, unlist_arandi, true_clusters[[d]])
    
    
    comparison_data_entry <- tibble(
      Dataset = d,
      # Clustering = list(.curr_clusters),
      Adjusted_rand_index = curr_rand_vec,
      Type = paste0(curr_nice_type) #,
      # Seed = k # ,
      # Weight = curr_weight
    )
    
    comparison_data <- rbind(comparison_data, comparison_data_entry)
  }
}

# === Plotting =================================================================

new_order <- unique(comparison_data$Type)[c(3, 1, 4, 2, 5, 7:14, 6)]
comparison_data_1 <- within(comparison_data, 
                            Type <- factor(Type, 
                                           levels=new_order)
)


bp <- ggplot(data = comparison_data_1, aes(x = Type, y = Adjusted_rand_index)) +
  geom_boxplot() +
  facet_wrap(~Dataset) +
  # theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  labs(
    title = "Simulation 1: comparison of consensus clusterings and individual chains to true clustering",
    y = "Adjusted rand index"
    
  )

plot_dir <- "Notes/Thesis/Images/Gen_data/Case_1"
dir.create(plot_dir, showWarnings = F)

ggsave(paste0(plot_dir, "/box_plot_ari_true_clustering_burn_in.png"), plot = bp)


# long_nice_types <- paste("Long run: seed", 1:10)

comparison_data_2 <- comparison_data
comparison_data_2$Type[comparison_data_2$Type %in% long_nice_types] <- "Long"
new_order_2 <- unique(comparison_data_2$Type)[c(3, 1, 4, 2, 5)]
comparison_data_2 <- within(comparison_data_2, 
                            Type <- factor(Type, 
                                           levels=new_order_2)
)


bp_collapsed <- ggplot(data = comparison_data_2, aes(x = Type, y = Adjusted_rand_index)) +
  geom_boxplot() +
  facet_wrap(~Dataset) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  labs(
    title = "Simulation 1: comparison of consensus clusterings and collapsed chains to true clustering",
    # subtitle = "Consensus clustering for various iterations, long runs corresponds to 10 chains of 2 million iterations",
    y = "Adjusted rand index"
    
  )
ggsave(paste0(plot_dir, "/box_plot_ari_true_clustering_collapsed_long_burn_in.png"), plot = bp_collapsed)



vp <- ggplot(data = comparison_data_2, aes(x = Type, y = Adjusted_rand_index)) +
  geom_violin() +
  facet_wrap(~Dataset) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  labs(
    title = "Simulation 1: comparison of consensus clusterings and collapsed chains to true clustering",
    # subtitle = "Consensus clustering for various iterations, long runs corresponds to 10 chains of 2 million iterations",
    y = "Adjusted rand index"
    
  )

ggsave(paste0(plot_dir, "/violin_plot_ari_true_clustering_collapsed_long_burn_in.png"), plot = vp)

