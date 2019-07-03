
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

# For controlling ordering of boxplots
add_level <- function(x, new_level){
  if(is.factor(x)) return(factor(x, levels=c(levels(x), new_level)))
  return(x)
}

# === Main script ==============================================================

# Set ggplot2 theme
theme_set(theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust=1)))

# The directroy holding the data
gen_dir <- "/home/MINTS/sdc56/Desktop/Gen_data_output/MDI_output"

# Directroy to save plots to
plot_dir <- "Notes/Thesis/Images/Gen_data/Case_2"
dir.create(plot_dir, showWarnings = F)

# The true labelling for this generated data
true_labels <- c(rep(1, 25),
                 rep(2, 50),
                 rep(3, 75),
                 rep(4, 100),
                 rep(5, 150)
)

# All directories holding MDI output tibbles
all_dirs_in_full <- list.dirs(path = gen_dir, recursive = F)
all_dirs_short <- list.dirs(path = gen_dir, recursive = F, full.names = F)

# Select the Long chain dirs
longer_dirs <- all_dirs_in_full %>% 
  extract(grep("Long" , .))

# select the consensus dirs
consensus_dirs <- all_dirs_in_full %>% 
  extract(grep("Consensus" , .)) 

tib_types <- list.dirs(path = gen_dir, recursive = F, full.names = F)

# gen_long_run_name <- c("Long", "Longer")
# gen_many_seeds_name <- "Many_seeds_"

# This is now redundant but goes back to bad folder organisation
all_rel_dirs <- all_dirs_in_full # c(consensus_dirs, longer_dirs)

# I don't think this is used now
types <- stringr::str_replace(tib_types, "_", " ")

# The name of the object the MDI output data is saved to in each dir in all_dirs_in_full
tib_name <- "compare_tibble.rds"
all_data <- paste0(all_rel_dirs, "/", tib_name)

# Read in the tibbles
all_tibs <- lapply(all_data, readRDS)

# Separate out the types of MDI applied
consensus_types <- tib_types %>% extract(grep("Consensus", .))
long_types <- tib_types %>% extract(grep("Long", .))

# Used in making nicer x-axis annotation
nice_types <- c(str_replace(consensus_types, "_", " "), str_replace(long_types, "_", " run: seed " ))


# === Set up to record =========================================================


# The object that will hold the comparison data
comparison_data <- tibble(
  Adjusted_rand_index = numeric(),
  Dataset = character(),
  Type = character()
)

# The object that will hold the current iteration of the loop
comparison_data_entry <- comparison_data

# Datasets present
datasets <- all_tibs[[1]]$dataset
n_datasets <- length(datasets)

# Number of MDI outputs here
n_mdi <- length(all_data)

# Loop over different MDI outputs and save comparison to true allocation
for(i in 1:n_mdi){
  curr_tib <- all_tibs[[i]]
  curr_type <- types[i]
  curr_nice_type <- nice_types[i]
for(d in datasets){
  curr_index <- which(curr_tib$dataset == d)
  
  curr_clustering <- curr_tib$mdi_allocation[[curr_index]]
  
  
  curr_rand_vec <- apply(curr_clustering, 1, unlist_arandi, true_labels)
  
  
  comparison_data_entry <- tibble(
    Adjusted_rand_index = curr_rand_vec,
    Dataset = d,
    Type = curr_nice_type
  )
  
  comparison_data <- rbind(comparison_data, comparison_data_entry)
}
}

# === Plotting =============================================================

# Find the appropriate order for boxplots
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
    title = "Generated data: comparison of consensus clusterings and individual chains",
    y = "Adjusted rand index"
    
  )


ggsave(paste0(plot_dir, "/box_plot_ari_true_clustering.png"), plot = bp)

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
    title = "Generated data: comparison of consensus clusterings and collapsed chains",
    # subtitle = "Consensus clustering for various iterations, long runs corresponds to 10 chains of 2 million iterations",
    y = "Adjusted rand index"
    
  )
ggsave(paste0(plot_dir, "/box_plot_ari_true_clustering_collapsed_long.png"), plot = bp_collapsed)


vp <- ggplot(data = comparison_data_2, aes(x = Type, y = Adjusted_rand_index)) +
  geom_violin() +
  facet_wrap(~Dataset) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  labs(
    title = "Generated data: comparison of consensus clusterings and collapsed chains",
    # subtitle = "Consensus clustering for various iterations, long runs corresponds to 10 chains of 2 million iterations",
    y = "Adjusted rand index"
    
  )

ggsave(paste0(plot_dir, "/violin_plot_ari_true_clustering_collapsed_long.png"), plot = vp)
