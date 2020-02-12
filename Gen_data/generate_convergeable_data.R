#!/usr/bin/env Rscript

# For tibbles
library(tibble) # for dataframe of lists

# For data wrangling
library(dplyr)

# Heatmapping
library(pheatmap) # install.packages("pheatmap", dep = T)

# Colour palettes
library(RColorBrewer)

# For command line arguments
library(optparse) # install.packages("optparse")

# Load data.table to access fread and fwrite functions
library(data.table) # install.packages("data.table", dep = T)

# Load magrittr for the pipe %>%
library(magrittr)

# for mvrnorm
library(MASS)

# Reproducibility
set.seed(1)

# Directory to save data to
save_dir <- "~/Documents/PhD/Year_1/Consensus_clustering/Data/Generated_data/Easier_convergence/"
file_name <- "Test_data_"
ph_name <- "heatmap_"
ph_title <- "Test dataset "

# Number of columns
d <- 40

# Number of clusters present
n_clust <- 11

# Means for each cluster
means <- seq(-floor(0.5 * n_clust), floor(0.5 * (n_clust - 1)))

# Number of points in each cluster
n <- rep(10, n_clust)

# Covariance matrix for all (in line with MDI assumptions)
my_cov <- diag(d)

# Row names
row_names <- paste0("Person_", 1:nrow(df))

# Object to hold cluster data
my_data <- list()

# Create list of cluster matrices
for(i in 1:n_clust){
  mu <- rep(means[i], d)
  my_data[[i]] <- mvrnorm(n[i], mu, my_cov)
}

# Combine these into one dataset
df <- do.call(rbind, my_data)
row.names(df) <- row_names

# Annotation data for pheatmap
# This data.frame describes how many points are associated with each cluster
label_df <- data.frame(Label = 1:n_clust, N = n)

# Vector to hold membership label
labels <- c()
for(i in 1:n_clust){
  labels <- c(labels, rep(label_df[i, 1], label_df[i, 2]))
}

# Put this into a data.frame for annotation purposes
annotation_labels <- data.frame(Labels = as.factor(labels))
row.names(annotation_labels) <- row_names

# Heatmap!
pheatmap(df, annotation_row = annotation_labels)

# Now create 3 sets of indices permuting the dataset so different points are in different clusters 
# To do this we permute row names
ind_1 <- 1:nrow(df)
ind_2 <- sample(1:nrow(df), nrow(df), replace = F)
ind_3 <- sample(ind_2, nrow(df), replace = F)

df_1 <- df %>% 
  set_rownames(row_names[ind_1])

annotation_labels_1 <- annotation_labels %>% 
  set_rownames(row_names[ind_1])

df_2 <- df %>% 
  set_rownames(row_names[ind_2]) %>% 
  .[match(row.names(df_1), row.names(df_2)), ]

annotation_labels_2 <- annotation_labels %>% 
  set_rownames(row_names[ind_2])

df_3 <- df %>% 
  set_rownames(row_names[ind_3]) %>% 
  .[match(row.names(df_1), row.names(df_3)), ]

annotation_labels_3 <- annotation_labels %>% 
  set_rownames(row_names[ind_3])

head(df_1)
head(df_2)
head(df_3)

# Visualise the data - all are the same bar row names
pheatmap(df_1,
         cluster_rows = F,
         annotation_row = annotation_labels_1)

pheatmap(df_2,
         cluster_rows = F,
         annotation_row = annotation_labels_2)

pheatmap(df_3,
         cluster_rows = F,
         annotation_row = annotation_labels_3)

pheatmap(df_1,
         cluster_rows = T,
         annotation_row = annotation_labels_1,
         main = paste0(ph_title, 1),
         filename = paste0(save_dir, ph_name, 1, ".png"))

pheatmap(df_2,
         cluster_rows = T,
         annotation_row = annotation_labels_2,
         main = paste0(ph_title, 2),
         filename = paste0(save_dir, ph_name, 2, ".png"))

pheatmap(df_3,
         cluster_rows = T,
         annotation_row = annotation_labels_3,
         main = paste0(ph_title, 3),
         filename = paste0(save_dir, ph_name, 3, ".png"))

write.csv(df_1, paste0(save_dir, file_name, "1.csv"))
write.csv(df_2, paste0(save_dir, file_name, "2.csv"))
write.csv(df_3, paste0(save_dir, file_name, "3.csv"))
