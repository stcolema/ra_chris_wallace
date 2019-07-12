
library(data.table)
library(pheatmap)
library(tibble)
library(magrittr)
library(ggplot2)
source("./Analysis/Analysis_script_functions/plot_comparison_expression_clustering.R")

# Set a seed for reproducibility
set.seed(1)


# Source the function to find gene sets from the KEGG data
source("./Analysis/read_in_kegg_sets.R")

# Read in the universe of genes for this project
probe_key <- fread("~/Desktop/Final_set/Meta_data/probe_key_unique_names.csv")
universe <- probe_key$Gene
my_universe <- probe_key$Unique_gene_name
reduced_universe <- unique(universe)

# The KEGG data
kegg_data <- "Data/kegg_msigdb.txt"

# The pathways to include
pathway_names <- c(
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_INOSITOL_PHOSPHATE_METABOLISM",
  "KEGG_INFLAMMATORY_BOWEL_DISEASE"
)

# This pathway is not present in the KEGG download so we had to download it 
# separately and parse it with python
missing_path <- pathway_names[3]

# Find the gene sets and reduce to those deisred
gene_sets <- find_gene_sets(kegg_data, universe) #, unique_universe = my_universe)
rel_gene_sets <- gene_sets[names(gene_sets) %in% pathway_names]

# The IBD set is missing from the above, read-in the file created by calling 
# python3 find_ibd_genes.py ./kegg_ibd_genes.txt ibd_genes.csv
ibd_set <- fread("Data/ibd_genes.csv", header = T) %>% unlist() %>% unname()

# Add to our list of gene sets using the unique labels
rel_gene_sets[[missing_path]] <- universe[universe %in% ibd_set]

# Reduce to a vector
genes_to_inclue <- rel_gene_sets %>% unlist() %>% unname() %>% unique()

genes_to_include_rep <- universe[universe %in% genes_to_inclue]
unique_genes_to_include <- my_universe[universe %in% genes_to_inclue]

other_genes_to_include <- universe[! universe %in% genes_to_inclue]
other_unique_genes_to_include <- my_universe[! universe %in% genes_to_inclue]

# Create a data.table that contains the validation information,
# i.e. the pathways the genes belong to
validation_dt <- data.table(Gene_ID = c(unique_genes_to_include, other_unique_genes_to_include),
                            Gene = c(genes_to_include_rep, other_genes_to_include),
                            KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY = FALSE,
                            KEGG_INOSITOL_PHOSPHATE_METABOLISM = FALSE,
                            KEGG_INFLAMMATORY_BOWEL_DISEASE = FALSE
                            )

for(path in pathway_names){
  path_genes <- rel_gene_sets[[path]]
  
  
  
  validation_dt[validation_dt[["Gene"]] %in% path_genes, path] <- T
}


# Location of tibble
# main_dir <- "~/Desktop/Final_set/CD/MDI_output/Consensus_1000/"
tibble_name <- "compare_tibble.rds"

# main_dir <- "~/Desktop/Output_324/Consensus_500_324/"
main_dir <- "~/Desktop/Newnet_254_output/Analysis/"

# Tibble containing similarity matrices
my_tib <- readRDS(paste0(main_dir, tibble_name))

curr_sim <- my_tib$similarity_matrix[[1]]
curr_corr <- my_tib$correlation_matrix[[1]]

probes_present_per_dataset <- read.csv("~/Desktop/newnet_228/Meta_data/probes_present_per_dataset.csv", row.names = 1, header = T)

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

new_row_names <- probe_key$Unique_gene_name[match(row.names(probes_present_per_dataset), probe_key$ProbeID)]
row.names(probes_present_per_dataset) <- new_row_names

annotation_row <- as.data.frame(validation_dt[, 3:5]) %>% 
  set_rownames(validation_dt$Gene_ID) %>% 
  extract(match(row.names(curr_sim), validation_dt$Gene_ID),)

str(annotation_row)
annotation_row %>% as.matrix()
as.matrix(annotation_row)

Var1 = c("navy", "darkgreen")
names(Var1) = c(T, F)
Var2 = c("lightgreen", "navy")
names(Var2) = c(T, F)
Var3 = c("blue", "orange")
names(Var3) = c(T, F)

ann_colors = list(KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY = Var1, 
                  KEGG_INOSITOL_PHOSPHATE_METABOLISM = Var2,
                  KEGG_INFLAMMATORY_BOWEL_DISEASE = Var3)

str(curr_sim)
str(annotation_row)


colnames(annotation_row) <- c(
  "NOD_LIKE",
  "INOSITOL",
  "IBD"
)

Var1 = c("blue", "orange")
names(Var1) = c("Member", "Non-Member")
Var2 = c("blue", "orange")
names(Var2) = c("Member", "Non-Member")
Var3 = c("blue", "orange")
names(Var3) = c("Member", "Non-Member")

Var2 <- Var1 <- Var3

ann_colors = list(NOD_LIKE = Var1, 
                  INOSITOL = Var2,
                  IBD = Var3)

annotation_row_2 <- annotation_row
annotation_row_2[annotation_row_2 == T] <- "Member"
annotation_row_2[annotation_row_2 == F] <- "Non-Member"

annotation_row_3 <- annotation_row_2[order(annotation_row_2$NOD_LIKE),]
curr_sim_2 <- curr_sim[order(annotation_row_2$NOD_LIKE),]

# annotation_row_4 <- annotation_row_2 %>%
#   inset2(., Rowname, row.names(annotation_row_2)) %>% 
#   dplyr::arrange(NOD_LIKE, INOSITOL)

annotation_row_4 <- annotation_row_2
annotation_row_4[["Rownames"]] <- row.names(annotation_row_4)
annotation_row_4 <- annotation_row_4 %>%   dplyr::arrange(NOD_LIKE, INOSITOL, IBD)
row.names(annotation_row_4) <- annotation_row_4$Rownames

annotation_row_4 <- annotation_row_4[,1:3]


curr_sim_3 <- curr_sim[match(row.names(annotation_row_4), row.names(curr_sim)), match(row.names(annotation_row_4), row.names(curr_sim))]

col_pal <- colorRampPalette(c("white", "#146EB4"))(100)
breaks_0_1 <- c(
  seq(0, 1, length.out = ceiling(length(col_pal)) + 1)
)

col_pal <- colorRampPalette(c("#FF9900", "white", "#146EB4"))(100)
cor_pal <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(100)
palette_length <- length(col_pal)

my_breaks <- c(
  seq(-1, 0, length.out = ceiling(palette_length / 2) + 1),
  seq(1 / palette_length, 1, length.out = floor(palette_length / 2))
)



cl_pred <- mcclust::maxpear(curr_sim)$cl
cl_true <- annotation_row_2$INOSITOL
cl_true[cl_true == "Member"] <- 1
cl_true[cl_true == "Non-Member"] <- 2
cl_true <- as.numeric(cl_true)
mcclust::arandi(cl_pred, cl_true)

# === Heatmapping ==============================================================

row_order <- hclust(dist(curr_sim))$order

ph1 <- pheatmap(curr_sim[row_order, row_order],
         annotation_row = annotation_row_2,
         annotation_colors = ann_colors,
         color = col_pal,
         breaks = my_breaks,
         show_rownames = F,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         silent = F)$gtable

# row_order <- ph1$tree_row$order

ph2 <- pheatmap(curr_corr[row_order, row_order],
                annotation_row = annotation_row_2,
                annotation_colors = ann_colors,
                color = cor_pal,
                breaks = my_breaks,
                show_rownames = F,
                show_colnames = F,
                cluster_rows = F,
                cluster_cols = F,
                silent = F)$gtable

ph_list <- list(ph1, ph2)

inostiol_genes <- rel_gene_sets$KEGG_INOSITOL_PHOSPHATE_METABOLISM

gene_ids <- probe_key$Unique_gene_name[probe_key$Gene %in% inostiol_genes]
indices_kept <- match(gene_ids, row.names(curr_sim))
sub_sim <- curr_sim[indices_kept, indices_kept]

sub_corr <- curr_corr[indices_kept, indices_kept]

sub_row_order <- hclust(dist(sub_sim))$order

ph_sub_1 <- pheatmap(sub_sim[sub_row_order, sub_row_order],
         annotation_row = annotation_row_2,
         annotation_colors = ann_colors,
         color = col_pal,
         breaks = my_breaks,
         show_rownames = F,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         silent = F)$gtable

ph_sub_2 <- pheatmap(sub_corr[sub_row_order, sub_row_order],
         annotation_row = annotation_row_2,
         annotation_colors = ann_colors,
         color = cor_pal,
         breaks = my_breaks,
         show_rownames = F,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         silent = F)$gtable


ph_sub_list <- list(ph_sub_1, ph_sub_2)

combine_pheatmaps(ph_list,
                  save_name = "~/Desktop/cd14_comp_psm_corr.png", 
                  main = "CD14: Comparison of annotated PSM and Correlation matrix")

combine_pheatmaps(ph_sub_list, 
                  save_name = "~/Desktop/cd14_comp_sub_psm_corr.png", 
                  main = "CD14: Comparison of annotated PSM and Correlation matrix (Inostiol)")

# === Investigate inostiol clustering probability ==============================

non_empty_inostiol_ind <- probes_present_per_dataset$CD14[match(row.names(sub_sim), row.names(probes_present_per_dataset))]
non_empty_inostiol_sim <- sub_sim[non_empty_inostiol_ind, non_empty_inostiol_ind]

inostiol_mean <- mean(non_empty_inostiol_sim[lower.tri(non_empty_inostiol_sim, diag = FALSE)])

find_mean_prob_distribution <- function(sim, n_sample, n_subset){
  
  n_genes <- nrow(sim)
  mean_prob <- c()
  for(i in 1:n_sample){
      ind_ <- sample(1:n_genes, replace = F, size = n_subset)
      sub_sim <- sim[ind_, ind_]
      mean_prob <- c(mean_prob, mean(sub_sim[lower.tri(sub_sim, diag = FALSE)]))
  }
  mean_prob
}

n_sample <- 1e5
m <- length(non_empty_inostiol_ind)

# non_empty_inostiol_ind <- probes_present_per_dataset$CD14[match(row.names(sub_sim), row.names(probes_present_per_dataset))]
no_inostiol_ind <- which(! row.names(curr_sim) %in% gene_ids)

sim_less_inostiol <- curr_sim[no_inostiol_ind, no_inostiol_ind]

non_empty_non_inostiol_ind <- probes_present_per_dataset$CD14[match(row.names(sim_less_inostiol), row.names(probes_present_per_dataset))]
non_empty_non_inostiol_sim <- sim_less_inostiol[non_empty_non_inostiol_ind, non_empty_non_inostiol_ind]

non_empty_mean_dstn <- find_mean_prob_distribution(non_empty_non_inostiol_sim, n_sample, m)
non_empty_non_inostiol_dt <- data.frame(Prob = non_empty_mean_dstn, Class = "Non-empty")

p_non_empty_non_inostiol <- ggplot(non_empty_non_inostiol_dt, aes(x = Prob, fill = Class)) +
  geom_density() +
  geom_vline(xintercept = inostiol_mean, linetype = "dashed", color = "black") +
  theme_bw() +
  labs(
    title = "Distribution of mean probability of pairwise alignment"
  )

p_non_empty_non_inostiol

# inostiol_mean <- mean(non_empty_non_inostiol_sim[lower.tri(non_empty_non_inostiol_sim, diag = FALSE)])

non_empty_ind <- which(probes_present_per_dataset$CD14[match(row.names(curr_sim), row.names(probes_present_per_dataset))])
non_empty_sim <- curr_sim[non_empty_ind, non_empty_ind]



mean_prob_non_empty <- find_mean_prob_distribution(non_empty_sim, n_sample, m)

my_probs <- data.frame(Prob = mean_prob_non_empty, Class = "Non-empty genes")

p_non_empty_all <- ggplot(my_probs, aes(x = Prob, fill = Class)) +
  geom_density() +
  geom_vline(xintercept = inostiol_mean, linetype = "dashed", color = "black") +
  theme_bw() +
  labs(
    title = "Distribution of mean probability of pairwise alignment"
  )

p_non_empty_all

# Empty included

# 
# inostiol_mean_incl_empty <- mean(sub_sim[lower.tri(sub_sim, diag = FALSE)])
# n_genes <- nrow(curr_sim)
# m <- nrow(sub_sim)
# mean_prob <- c()
# for(i in 1:10000){
#   ind_ <- sample(1:n_genes, replace = F, size = m)
#   curr_sub_sim <- curr_sim[ind_, ind_]
#   mean_prob <- c(mean_prob, mean(curr_sub_sim[lower.tri(curr_sub_sim, diag = FALSE)]))
# }
# 
# sum(inostiol_mean_incl_empty > mean_prob) / length(mean_prob)
# 
# plot(density(mean_prob)) 
# 
# my_probs_2 <- data.frame(Prob = mean_prob, Class = "Including absent genes")
# 
# my_probs_df <- rbind(my_probs, my_probs_2)
# 
# ggplot(my_probs_df, aes(x = Prob, fill = Class, alpha =0.4)) +
#   geom_density() +
#   geom_vline(xintercept = inostiol_mean, linetype = "dashed", color = "blue") +
#   geom_vline(xintercept = inostiol_mean_incl_empty, linetype = "dashed", color = "red") +
#   theme_bw() +
#   labs(
#     title = "Distribution of mean probability of pairwise alignment"
#   )

# === Distribution of PSM values ===============================================

plot(density(non_empty_sim[lower.tri(non_empty_sim, diag = F)]))
plot(density(non_empty_non_inostiol_sim[lower.tri(non_empty_non_inostiol_sim, diag = F)]))
plot(density(non_empty_inostiol_sim[lower.tri(non_empty_inostiol_sim, diag = F)]))

non_inostiol_psm_data_frame <- data.frame(Proportion = non_empty_non_inostiol_sim[lower.tri(non_empty_non_inostiol_sim, diag = F)],
Class = "Other")

inostiol_psm_data_frame <-data.frame(Proportion = non_empty_inostiol_sim[lower.tri(non_empty_inostiol_sim, diag = F)],
                                         Class = "Inostiol") 

non_empty_psm_data_frame <- rbind(inostiol_psm_data_frame, non_inostiol_psm_data_frame)


violin_all <- ggplot(non_empty_psm_data_frame, aes(x = Class, y = Proportion, fill = Class)) +
  geom_violin() +
  theme_bw() +
  labs(
    title = "Density plots for PSM entries of Inostiol and other genes"
  )

boxplot_all <- ggplot(non_empty_psm_data_frame, aes(x = Class, y = Proportion, fill = Class)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = "Density plots for PSM entries of Inostiol and other genes"
  )

ggsave("~/Desktop/psm_violin.png", plot = violin_all)
ggsave("~/Desktop/psm_boxplot.png", plot = boxplot_all)

all_psm_density <- ggplot(non_empty_psm_data_frame, aes(x = Proportion, fill = Class, alpha =0.4)) +
  geom_density() +
  theme_bw() +
  labs(
    title = "Density plots for PSM entries of Inostiol and other genes"
  )

ggsave("~/Desktop/psm_density.png", plot = all_psm_density)

inostiol_psm_density <- ggplot(inostiol_psm_data_frame, aes(x = Proportion, fill = Class)) +
  geom_density() +
  theme_bw() +
  labs(
    title = "Density plots for PSM entries of Inostiol genes"
  )

ggsave("~/Desktop/inotiol_psm_density.png", plot = inostiol_psm_density)

non_inostiol_psm_density <- ggplot(non_inostiol_psm_data_frame, aes(x = Proportion, fill = Class)) +
  geom_density() +
  theme_bw() +
  labs(
    title = "Density plots for PSM entries of other genes"
  ) +
  scale_fill_manual(values=c("#999999"))

ggsave("~/Desktop/non_inotiol_psm_density.png", plot = non_inostiol_psm_density)

  

# === Subsection of inostiol genes in correlation and PSM ======================

sub_sub_sim <- sub_sim[sub_row_order, sub_row_order]
sub_sub_corr <- sub_corr[sub_row_order, sub_row_order]

ph_sub_sub_1 <- pheatmap(sub_sub_sim[-c(37:60),-c(37:60)],
                     annotation_row = annotation_row_2,
                     annotation_colors = ann_colors,
                     color = col_pal,
                     breaks = my_breaks,
                     show_rownames = F,
                     show_colnames = F,
                     cluster_rows = F,
                     cluster_cols = F,
                     silent = F)$gtable




ph_sub_sub_2 <- pheatmap(sub_sub_corr[-c(37:60),-c(37:60)],
                     annotation_row = annotation_row_2,
                     annotation_colors = ann_colors,
                     color = cor_pal,
                     breaks = my_breaks,
                     show_rownames = F,
                     show_colnames = F,
                     cluster_rows = F,
                     cluster_cols = F,
                     silent = F)$gtable

ph_sub_sub_list <- list(ph_sub_sub_1, ph_sub_sub_2)

combine_pheatmaps(ph_sub_sub_list, 
                  save_name = "~/Desktop/cd14_comp_sub_sub_psm_corr.png", 
                  main = "CD14: PSM and Correlation matrix (Inostiol, non-empty)")


sub_cor_order <- hclust(dist(sub_sub_corr[-c(37:60),-c(37:60)]))$order

non_empty_inostiol_sim <- sub_sub_sim[-c(37:60),-c(37:60)]
non_empty_inostiol_corr <- sub_sub_corr[-c(37:60),-c(37:60)]

ph_sub_sub_1_corr_ord <- pheatmap(non_empty_inostiol_sim[sub_cor_order, sub_cor_order],
                         annotation_row = annotation_row_2,
                         annotation_colors = ann_colors,
                         color = col_pal,
                         breaks = my_breaks,
                         show_rownames = F,
                         show_colnames = F,
                         cluster_rows = F,
                         cluster_cols = F,
                         silent = F)$gtable

ph_sub_sub_2_corr_ord <- pheatmap(non_empty_inostiol_corr[sub_cor_order, sub_cor_order],
                         annotation_row = annotation_row_2,
                         annotation_colors = ann_colors,
                         color = cor_pal,
                         breaks = my_breaks,
                         show_rownames = F,
                         show_colnames = F,
                         cluster_rows = F,
                         cluster_cols = F,
                         silent = F)$gtable

ph_sub_sub_list_corr_ord <- list(ph_sub_sub_1_corr_ord, ph_sub_sub_2_corr_ord)

combine_pheatmaps(ph_sub_sub_list_corr_ord, 
                  save_name = "~/Desktop/cd14_comp_sub_sub_psm_corr_order.png", 
                  main = "CD14: PSM and Correlation matrix (Inostiol, non-empty)")


pheatmap(curr_sim_2,
         annotation_row = annotation_row_3,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         color = col_pal,
         breaks = my_breaks)

pheatmap(curr_sim_3,
         annotation_row = annotation_row_4,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         color = col_pal,
         breaks = my_breaks)

ugly <- cbind(annotation_row_4, curr_sim_3)

sub_1_curr_sim_3 <- curr_sim_3[annotation_row_4$NOD_LIKE == "Member" & annotation_row_4$IBD == "Member", ]

hc_1 <- hclust(dist(sub_1_curr_sim_3))
order_1 <- hc_1$order

n_1 <- nrow(sub_1_curr_sim_3)

sub_2_curr_sim_3 <- curr_sim_3[annotation_row_4$NOD_LIKE == "Member" & annotation_row_4$IBD == "Non-Member", ]

hc_2 <- hclust(dist(sub_2_curr_sim_3))
order_2 <- hc_2$order + n_1

n_2 <- nrow(sub_2_curr_sim_3) + n_1

sub_3_curr_sim_3 <- curr_sim_3[annotation_row_4$INOSITOL == "Member", ]

hc_3 <- hclust(dist(sub_3_curr_sim_3))
order_3 <- hc_3$order + n_2

n_3 <- nrow(sub_3_curr_sim_3) + n_2

sub_4_curr_sim_3 <- curr_sim_3[annotation_row_4$NOD_LIKE == "Non-Member" & annotation_row_4$IBD == "Member", ]

hc_4 <- hclust(dist(sub_4_curr_sim_3))
order_4 <- hc_4$order + n_3

n_4 <- nrow(sub_4_curr_sim_3) + n_3

sub_5_curr_sim_3 <- curr_sim_3[annotation_row_4$NOD_LIKE == "Non-Member" & annotation_row_4$IBD == "Non-Member" & annotation_row_4$INOSITOL == "Non-Member", ]

hc_5 <- hclust(dist(sub_5_curr_sim_3))
order_5 <- hc_5$order + n_4

n_5 <- nrow(sub_5_curr_sim_3) + n_4

new_order <- c(order_1, order_2, order_3, order_4, order_5)

pheatmap(curr_sim_3[new_order, new_order],
         annotation_row = annotation_row_4[new_order,],
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         color = col_pal,
         breaks = breaks_0_1)

pheatmap(ugly)

# my_ann <- data.frame(1:1000 %% 2) %>% set_rownames(row.names(curr_sim))

# pheatmap(curr_sim,
#          annotation_row = my_ann)
# ,
         # annotation_colors = ann_colors)
