#!/usr/bin/env Rscript

# === Libraries ================================================================

library(data.table)
library(pheatmap)
library(tibble)
library(magrittr)
library(ggplot2)
library(ggforce)
library(plyr)
library(optparse)
library(rlist)

source("~/Desktop/ra_chris_wallace/Analysis/Analysis_script_functions/plot_comparison_expression_clustering.R")

# Source the function to find gene sets from the KEGG data
source("~/Desktop/ra_chris_wallace//Analysis/read_in_kegg_sets.R")


# Set a seed for reproducibility
set.seed(1)

# === Functions ================================================================

# Returns a list of the gtable of pheatmaps with common order for sim and cor_mat
sim_cor_ph_list <- function(sim, cor_mat, annotation_row, annotation_colours,
                            col_pal = colorRampPalette(c("#FF9900", "white", "#146EB4"))(100),
                            cor_pal = colorRampPalette(c("#146EB4", "white", "#FF9900"))(100),
                            breaks = NULL) {

  # Set up breaks if not defined
  if (is.null(breaks)) {
    palette_length <- length(col_pal)
    breaks <- c(
      seq(-1, 0, length.out = ceiling(palette_length / 2) + 1),
      seq(1 / palette_length, 1, length.out = floor(palette_length / 2))
    )
  }

  # Order by the PSM
  row_order <- hclust(dist(sim))$order

  # Compare heatmaps of PSM and Correlation matrix for entire dataset
  ph1 <- pheatmap(sim[row_order, row_order],
    annotation_row = annotation_row,
    annotation_colors = annotation_colours,
    color = col_pal,
    breaks = my_breaks,
    show_rownames = F,
    show_colnames = F,
    cluster_rows = F,
    cluster_cols = F,
    silent = T
  )$gtable

  show_rownames <- T
  if (nrow(sim) > 70) {
    show_rownames <- F
  }

  ph2 <- pheatmap(cor_mat[row_order, row_order],
    annotation_row = annotation_row,
    annotation_colors = annotation_colours,
    color = cor_pal,
    breaks = my_breaks,
    show_rownames = show_rownames,
    show_colnames = F,
    cluster_rows = F,
    cluster_cols = F,
    silent = T
  )$gtable

  ph_list <- list(ph1, ph2)
}


# Find the distribution of the means of the probability of clsutering
# together for sets of n_subset genes
find_mean_prob_distribution <- function(sim, n_sample, n_subset) {
  n_genes <- nrow(sim)
  mean_prob <- c()
  for (i in 1:n_sample) {
    ind_ <- sample(1:n_genes, replace = F, size = n_subset)
    sub_sim <- sim[ind_, ind_]
    new_mean <- mean(sub_sim[lower.tri(sub_sim, diag = FALSE)])
    mean_prob <- c(mean_prob, new_mean)
  }
  mean_prob
}

# For a given PSM find the gene sets defined by the following:
#   1. Find genes with probability of being assigned to the same cluster above a
# threshold
#   2. Check if either of the genes is already in a set; if so assign them to this set
#   3. Optionally combine sets
find_gene_combinations <- function(sim, threshold = 0.8) { # , combine = TRUE) {
  indices_to_inspect_all <- which(sim > threshold, arr.ind = TRUE)

  indices_to_inspect <- indices_to_inspect_all[indices_to_inspect_all[, 1] != indices_to_inspect_all[, 2], ]

  p_combination <- list()

  # res <- lapply(p_combination, function(ch) grep("PLCB2", ch))
  # sapply(res, function(x) length(x) > 0)
  count <- 0
  for (i in 1:nrow(indices_to_inspect)) {
    entries_to_drop <- c()

    curr_col <- indices_to_inspect[i, 2]
    curr_row <- indices_to_inspect[i, 1]
    p1 <- row.names(sim)[curr_row]
    p2 <- colnames(sim)[curr_col]

    cluster_to_join_1 <- lapply(p_combination, function(ch) grep(p1, ch))
    cluster_to_join_2 <- lapply(p_combination, function(ch) grep(p2, ch))

    ind_1_bool <- sapply(cluster_to_join_1, function(x) length(x) > 0)
    ind_2_bool <- sapply(cluster_to_join_2, function(x) length(x) > 0)

    if (length(ind_1_bool) > 0) {
      ind_1 <- which(ind_1_bool)
    }
    if (length(ind_2_bool) > 0) {
      ind_2 <- which(ind_2_bool)
    }

    if (length(ind_2_bool) > 0 & length(ind_1_bool) > 0) {
      no_cluster_existing <- sum(ind_1_bool, ind_2_bool) == 0
    } else {
      no_cluster_existing <- T
    }
    if (no_cluster_existing) {
      count <- count + 1
      p_combination[[count]] <- c(p1, p2)
    } else {
      if (length(ind_1) > 0) {
        for (j in ind_1) {
          p_combination[[j]] <- unique(c(p1, p2, p_combination[[j]]))
          if (length(ind_1) > 1) {
            p_combination[[ind_1[1]]] <- unique(c(p_combination[[ind_1[1]]], p_combination[[j]]))
            entries_to_drop <- c(entries_to_drop, j)
          }
        }
      }
      if (length(ind_2) > 0) {
        for (j in ind_2) {
          p_combination[[j]] <- unique(c(p1, p2, p_combination[[j]]))

          if (length(ind_2) > 1) {
            p_combination[[ind_2[1]]] <- unique(c(p_combination[[ind_2[1]]], p_combination[[j]]))
            entries_to_drop <- c(entries_to_drop, j)
          }
        }
      }
      entries_to_drop <- unique(entries_to_drop)

      for (j in entries_to_drop) {
        p_combination <- list.remove(p_combination, range = j)
      }
    }
  }
  compact(p_combination)
}




input_arguments <- function() {
  option_list <- list(

    # Directory to read from
    optparse::make_option(c("-d", "--dir"),
      type = "character",
      default = ".",
      help = "directory to read files from [default= %default]",
      metavar = "character"
    ),

    optparse::make_option(c("-p", "--probes"),
      type = "character",
      default = "~/Desktop/newnet_228/Meta_data/probes_present_per_dataset.csv",
      help = "file to read probes preesence across datasets from [default= %default]",
      metavar = "character"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

# === Set up inputs ============================================================

# Location of tibble
# main_dir <- "~/Desktop/Final_set/CD/MDI_output/Consensus_1000/"
tibble_name <- "compare_tibble.rds"

# main_dir <- "~/Desktop/Output_324/Consensus_500_324/"
# main_dir <- save_dir <- "~/Desktop/Newnet_254_output/Analysis/"
args <- input_arguments()
main_dir <- save_dir <- args$dir

# Plot type
plot_type <- ".png"

# Tibble containing similarity matrices
my_tib <- readRDS(paste0(main_dir, tibble_name))

# Palettes and breaks for heatmaps
col_pal <- colorRampPalette(c("#FF9900", "white", "#146EB4"))(100)
cor_pal <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(100)

# Read in the universe of genes for this project
probe_key <- fread("~/Desktop/Final_set/Meta_data/probe_key_unique_names.csv")
universe <- probe_key$Gene
my_universe <- probe_key$Unique_gene_name
reduced_universe <- unique(universe)

# The KEGG data
kegg_data <- "~/Desktop/ra_chris_wallace/Data/kegg_msigdb.txt"

# The pathways to include
pathway_names <- c(
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_INOSITOL_PHOSPHATE_METABOLISM",
  "KEGG_INFLAMMATORY_BOWEL_DISEASE"
)

abbreviated_names <- c(
  "NOD like",
  "Inostiol",
  "IBD"
)

n_pathways <- length(pathway_names)

# Probes present across datasets
probes_present_per_dataset_file_name <- args$probes
probes_present_per_dataset <- read.csv(probes_present_per_dataset_file_name,
  row.names = 1,
  header = T
)

colnames(probes_present_per_dataset) <- c(
  "CD14",
  "CD15",
  "CD19",
  "CD4",
  "CD8",
  "IL",
  "PLA",
  "RE",
  "TR"
)

# === Initialisation ===========================================================

palette_length <- length(col_pal)
my_breaks <- c(
  seq(-1, 0, length.out = ceiling(palette_length / 2) + 1),
  seq(1 / palette_length, 1, length.out = floor(palette_length / 2))
)

# Directory for heatmaps
save_dir_ph <- paste0(save_dir, "Heatmaps/")
dir.create(save_dir_ph, showWarnings = F)

# Subdirectories for heatmaps
save_dir_ph_all <- paste0(save_dir_ph, "All/")
dir.create(save_dir_ph_all, showWarnings = F)

save_dir_pathways <- vector("list", n_pathways) %>%
  set_names(pathway_names)

for (p in pathway_names) {
  save_dir_pathways[[p]] <- curr_dir <- paste0(save_dir_ph, p, "/")
  dir.create(curr_dir, showWarnings = F)
}

# === Update meta data =========================================================

# Update the row names of the probe key to genes
new_row_names <- probe_key$Unique_gene_name[match(row.names(probes_present_per_dataset), probe_key$ProbeID)]
row.names(probes_present_per_dataset) <- new_row_names

# === KEGG pathways ============================================================

# This pathway is not present in the KEGG download so we had to download it
# separately and parse it with python
missing_path <- pathway_names[3]

# Find the gene sets and reduce to those deisred
gene_sets <- find_gene_sets(kegg_data, universe) # , unique_universe = my_universe)
rel_gene_sets <- gene_sets[names(gene_sets) %in% pathway_names]

# The IBD set is missing from the above, read-in the file created by calling
# python3 find_ibd_genes.py ./kegg_ibd_genes.txt ibd_genes.csv
ibd_set <- fread("~/Desktop/ra_chris_wallace/Data/ibd_genes.csv", header = T) %>%
  unlist() %>%
  unname()

# Add to our list of gene sets using the unique labels
rel_gene_sets[[missing_path]] <- universe[universe %in% ibd_set]

# Reduce to a vector
genes_to_inclue <- rel_gene_sets %>%
  unlist() %>%
  unname() %>%
  unique()

# Find the gene ids of interest
genes_to_include_rep <- universe[universe %in% genes_to_inclue]
unique_genes_to_include <- my_universe[universe %in% genes_to_inclue]

other_genes_to_include <- universe[!universe %in% genes_to_inclue]
other_unique_genes_to_include <- my_universe[!universe %in% genes_to_inclue]


# Create a data.table that contains the validation information,
# i.e. the pathways the genes belong to
validation_dt <- data.table(
  Gene_ID = c(unique_genes_to_include, other_unique_genes_to_include),
  Gene = c(genes_to_include_rep, other_genes_to_include),
  KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY = FALSE,
  KEGG_INOSITOL_PHOSPHATE_METABOLISM = FALSE,
  KEGG_INFLAMMATORY_BOWEL_DISEASE = FALSE
)

for (path in pathway_names) {
  path_genes <- rel_gene_sets[[path]]

  validation_dt[validation_dt[["Gene"]] %in% path_genes, path] <- T
}

# === Preparing for heatmapping ================================================

# Current posterior similarity matrix
curr_sim <- my_tib$similarity_matrix[[1]]

# Create an annotation dataframe based on gene sets for heatmapping
annotation_row_logical <- as.data.frame(validation_dt[, 3:5]) %>%
  set_rownames(validation_dt$Gene_ID) %>%
  extract(match(row.names(curr_sim), validation_dt$Gene_ID), )

# Use strings (required for heatmap annotation)
annotation_row <- annotation_row_logical
annotation_row[annotation_row == T] <- "Member"
annotation_row[annotation_row == F] <- "Non-Member"

# Annotation colours indicating membership in a pathway

# Ideally we would use full pathway names, but these are too much for pheatmap
Var1 <- c("darkgreen", "white")
names(Var1) <- c("Member", "Non-Member")

# ann_colors = list(KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY = Var1,
#                   KEGG_INOSITOL_PHOSPHATE_METABOLISM = Var1,
#                   KEGG_INFLAMMATORY_BOWEL_DISEASE = Var1)

# Use abbreviated pathway names
colnames(annotation_row) <- c(
  abbreviated_names
)

ann_colors <- list(Var1, Var1, Var1) %>% set_names(abbreviated_names)

# === Data =====================================================================

# Datasets present
datasets <- my_tib$dataset
n_datasets <- length(datasets)

# Variables for FOR-loop
# d <- datasets[1]

# Current entry in tibble
# curr_ind <- which(my_tib$dataset == d)

# # Current posterior similarity matrix
# curr_sim <- my_tib$similarity_matrix[[1]]

# Current correlation matrix
# curr_corr <- my_tib$correlation_matrix[[curr_ind]]

# Pathway specific genes
pathway_gene_ids <- vector("list", n_pathways) %>%
  set_names(pathway_names)
for (p in pathway_names) {
  pathway_genes <- rel_gene_sets[[p]]
  pathway_gene_ids[[p]] <- probe_key$Unique_gene_name[probe_key$Gene %in% pathway_genes]
}

# === Heatmapping ==============================================================

# For each dataset create a comparison of annotated heatmaps for the PSM and
# correlation matrix for the entire dataset and then for each of the pathway
# specific genes

for (d in datasets) {

  # Current entry in tibble
  curr_ind <- which(my_tib$dataset == d)

  # Current posterior similarity matrix
  curr_sim <- my_tib$similarity_matrix[[curr_ind]]

  # Current correlation matrix
  curr_corr <- my_tib$correlation_matrix[[curr_ind]]

  ph_list <- sim_cor_ph_list(curr_sim, curr_corr, annotation_row, ann_colors,
    col_pal = col_pal,
    cor_pal = cor_pal,
    breaks = my_breaks
  )

  # Save name and main title
  save_all <- paste0(save_dir_ph_all, d, "_comp_psm_corr", plot_type)
  main_all <- paste0(d, ": Comparison of annotated PSM and Correlation matrix")

  # Combine and save heatmaps
  combine_pheatmaps(ph_list,
    save_name = save_all,
    main = main_all
  )

  # Pathway specific analysis
  for (i in 1:n_pathways) {
    p <- pathway_names[[i]]
    abrv_p <- abbreviated_names[[i]]

    # Current save directory
    curr_save_dir <- save_dir_pathways[[p]]

    # Current pathway gene ids
    curr_gene_ids <- pathway_gene_ids[[p]]

    # Subset similarity matrix and correlation matrix
    indices_kept <- match(curr_gene_ids, row.names(curr_sim))
    sub_sim <- curr_sim[indices_kept, indices_kept]
    sub_corr <- curr_corr[indices_kept, indices_kept]

    # Pheatmaps
    ph_sub_list <- sim_cor_ph_list(sub_sim, sub_corr, annotation_row, ann_colors,
      col_pal = col_pal,
      cor_pal = cor_pal,
      breaks = my_breaks
    )

    # Save name for heatmaps
    save_pathway <- paste0(curr_save_dir, d, "_comp_psm_corr", plot_type)

    # Title for heatmaps
    main_pathway <- paste0(d, ": Comparison of annotated PSM and Correlation matrix")

    # Combine and save heatmaps
    combine_pheatmaps(ph_sub_list,
      save_name = save_pathway,
      main = main_pathway
    )
  }
}

# === Investigate inostiol clustering probability ==============================

big_non_empty_non_pathway_dt <- data.frame(
  Prob = numeric(),
  Class = character(),
  Dataset = character(),
  Pathway = character(),
  Pathway_mean = numeric()
)

pathway_plots <- vector("list", n_datasets) %>% set_names(datasets)
n_sample <- 1e4

align_plot_dir <- paste0(save_dir, "Mean_alignment_probability/")
dir.create(align_plot_dir, showWarnings = F)

for (i in 1:n_datasets) {
  d <- datasets[i]

  align_plot_dataset_dir <- paste0(align_plot_dir) # , d, "/")
  dir.create(align_plot_dataset_dir, showWarnings = F)


  pathway_plots[[d]] <- vector("list", n_pathways) %>%
    set_names(pathway_names)

  # PSM for the current dataset
  curr_sim <- my_tib$similarity_matrix[[i]]

  # Pathway specific analysis
  for (j in 1:n_pathways) {
    p <- pathway_names[[j]]
    abrv_p <- abbreviated_names[[j]]
    # p <- pathway_names[[1]]

    # Current pathway gene ids
    curr_gene_ids <- pathway_gene_ids[[p]]

    # Subset similarity matrix and correlation matrix
    indices_kept <- match(curr_gene_ids, row.names(curr_sim))
    sub_sim <- curr_sim[indices_kept, indices_kept]

    # Find the non-empty probes / genes
    curr_probes_present <- probes_present_per_dataset[[d]]

    # Find the info about probe presence specific to the pathway
    non_empty_pathway_ind <- curr_probes_present[match(row.names(sub_sim), row.names(probes_present_per_dataset))]

    # Subset the PSM
    non_empty_pathway_sim <- sub_sim[non_empty_pathway_ind, non_empty_pathway_ind]

    # The mean probability of any pair of genes from the pathway clustering together
    pathway_mean <- mean(non_empty_pathway_sim[lower.tri(non_empty_pathway_sim, diag = FALSE)])

    # Number of samples to use to describe the distribution of means
    m <- sum(non_empty_pathway_ind)

    # Find the genes present in the PSM that are not members of the pathway set
    not_curr_path_ind <- which(!row.names(curr_sim) %in% curr_gene_ids)
    sim_less_curr_path <- curr_sim[not_curr_path_ind, not_curr_path_ind]

    # Find the non-empty genes not in the current pathway
    non_empty_non_pathway_ind <- curr_probes_present[match(row.names(sim_less_curr_path), row.names(probes_present_per_dataset))]
    non_empty_non_pathway_sim <- sim_less_curr_path[non_empty_non_pathway_ind, non_empty_non_pathway_ind]

    # Check if the number of non-empty non-pathway genes (m_2) is less than the
    # number of non-empty pathway genes (m). If this is the case we cannot
    # sample m genes from a set of m_2 genes
    m_2 <- nrow(non_empty_non_pathway_sim)
    if (m > m_2) {
      cat("Skipping", p, "for", d, "\n")
      next
    }

    # print(head(non_empty_non_pathway_sim))

    non_empty_mean_dstn <- find_mean_prob_distribution(non_empty_non_pathway_sim, n_sample, m)

    non_empty_non_pathway_dt <- data.frame(
      Prob = non_empty_mean_dstn,
      Class = "Non-empty",
      Dataset = d,
      Pathway = abrv_p,
      Pathway_mean = pathway_mean
    )

    # print(head(non_empty_non_pathway_dt))

    plot_title <- paste0(d, ": Distribution of mean probability of pairwise alignment (", abrv_p, ")")

    p_non_empty_non_pathway <- ggplot(non_empty_non_pathway_dt, aes(x = Prob, fill = Class)) +
      geom_density() +
      geom_vline(xintercept = pathway_mean, linetype = "dashed", color = "black") +
      theme_bw() +
      labs(
        title = plot_title
      )

    pathway_plots[[d]][[p]] <- p_non_empty_non_pathway

    # p_non_empty_non_pathway

    save_name <- paste0(align_plot_dataset_dir, d, "_", p, plot_type)

    ggsave(save_name, p_non_empty_non_pathway)

    big_non_empty_non_pathway_dt <- rbind(big_non_empty_non_pathway_dt, non_empty_non_pathway_dt)
  }
}

# === PSM violin plots =========================================================

density_save_dir <- paste0(save_dir, "PSM_densities/")
dir.create(density_save_dir, showWarnings = F)

gene_combinations <- vector("list", n_datasets) %>% set_names(datasets)

for (i in 1:n_datasets) {
  d <- datasets[i]
  curr_sim <- my_tib$similarity_matrix[[i]]

  gene_combinations[[d]] <- vector("list", n_pathways) %>% set_names(pathway_names)
  for (j in 1:n_pathways) {
    p <- pathway_names[[j]]
    abrv_p <- abbreviated_names[[j]]

    pathway_density_save_dir <- paste0(density_save_dir, p, "/")
    dir.create(pathway_density_save_dir, showWarnings = F)

    # p <- pathway_names[[1]]

    # Current pathway gene ids
    curr_gene_ids <- pathway_gene_ids[[p]]

    # Subset similarity matrix and correlation matrix
    indices_kept <- match(curr_gene_ids, row.names(curr_sim))
    sub_sim <- curr_sim[indices_kept, indices_kept]

    # Find the non-empty probes / genes
    curr_probes_present <- probes_present_per_dataset[[d]]

    # Find the info about probe presence specific to the pathway
    non_empty_pathway_ind <- curr_probes_present[match(row.names(sub_sim), row.names(probes_present_per_dataset))]

    # Subset the PSM
    non_empty_pathway_sim <- sub_sim[non_empty_pathway_ind, non_empty_pathway_ind]

    # Find the genes present in the PSM that are not members of the pathway set
    not_curr_path_ind <- which(!row.names(curr_sim) %in% curr_gene_ids)
    sim_less_curr_path <- curr_sim[not_curr_path_ind, not_curr_path_ind]

    # Find the non-empty genes not in the current pathway
    non_empty_non_pathway_ind <- curr_probes_present[match(row.names(sim_less_curr_path), row.names(probes_present_per_dataset))]
    non_empty_non_pathway_sim <- sim_less_curr_path[non_empty_non_pathway_ind, non_empty_non_pathway_ind]


    non_pathway_psm_data_frame <- data.frame(
      Proportion = non_empty_non_pathway_sim[lower.tri(non_empty_non_pathway_sim, diag = F)],
      Class = "Other"
    )

    pathway_psm_data_frame <- data.frame(
      Proportion = non_empty_pathway_sim[lower.tri(non_empty_pathway_sim, diag = F)],
      Class = abrv_p
    )

    non_empty_psm_data_frame <- rbind(pathway_psm_data_frame, non_pathway_psm_data_frame)

    plot_title <- paste0(d, ": Violin plots comparing PSM entries of ", abrv_p, " genes and all other genes")
    save_name <- paste0(pathway_density_save_dir, d, plot_type)

    violin_all <- ggplot(non_empty_psm_data_frame, aes(x = Class, y = Proportion, fill = Class)) +
      geom_violin() +
      theme_bw() +
      labs(
        title = plot_title
      )
    ggsave(save_name, violin_all)

    gene_combinations[[d]][[p]] <- find_gene_combinations(non_empty_pathway_sim, threshold = 0.8)
  }
}

gene_comb_save_dir <- paste0(save_dir, "/Gene_combinations/")
dir.create(gene_comb_save_dir, showWarnings = F)

for (d in datasets) {
  for (j in 1:n_pathways) {
    p <- pathway_names[j]
    save_name <- paste0(gene_comb_save_dir, d, "_", p, "_gene_combinations.txt")
    lapply(gene_combinations[[d]][[p]],
      write.table,
      save_name,
      append = T, sep = ", ", row.names = F
    )
  }
}

# lapply(gene_combinations$CD14$KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY,
#        write.table,
#        "~/Desktop/my_file.txt",
#        append = T, sep = ", ", row.names = F)
