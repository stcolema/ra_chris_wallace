#!/usr/env/bin/Rscript

################################################################################
#                                                                              #
# Title: Consensus inference performance                                       #
#                                                                              #
# Aim: Analyse saved clustering output of Mason's implementation of MDI for    #
# Consensus inference simulations run April 2020.                              #
# Compare models to the ground truth using ARI and a normalised Frobenius      #
# product.                                                                     #
#                                                                              #
# Author: Stephen Coleman                                                      #
# Date: 21/04/2020                                                             #
#                                                                              #
################################################################################

# For my ggplot theme settings
library(mdiHelpR)

# File name interactions
library(stringr)

# Plotting
library(ggplot2)

# I don't think this is used
library(tibble)

# To read in data
library(data.table)

# Visualising heatmaps (no longer used)
library(pheatmap)

# For frobenius.prod
library(matrixcalc)

# For arandi and maxpear for model evaluation and predicted clustering
library(mcclust)

# For facet_wrap_paginate
library(ggforce)

# For the saving of Sparse matrices
library(Matrix)

# For command line arguments
library(optparse)

# For pipe and related functions
library(magrittr)

my_dir <- "/Users/stephen/Desktop/simple_2d/Consensus/"

sim_results_files <- list.files(my_dir, pattern = ".csv", full.names = T)
n_files <- length(sim_results_files)

for(i in 1:n_files){
  x <- fread(sim_results_files[i], drop=1)
  
  if(i == 1){
    all_results_df <- x
  } else {
    all_results_df <- rbind(all_results_df, x)
  }
}

all_results_df$Simulation <- factor(all_results_df$Simulation)

# Set labels for facet wrapping
iter_labels <- c(paste0("Number of iterations: ", results_df$N_iter))
names(iter_labels) <- results_df$N_iter

seed_labels <- c(paste0("Number of chains: ", results_df$N_seeds))
names(seed_labels) <- results_df$N_seeds

# Plots
p1 <- all_results_df %>%
  ggplot(aes(x = N_seeds, y = ARI)) +
  geom_line(aes(group = Simulation), colour = "grey", alpha = 0.3) +
  facet_wrap(~N_iter, labeller = labeller(N_iter = iter_labels), ncol = 1) +
geom_smooth(color = "blue", lty = 1,
stat = "summary",
fill = "red",
alpha = 0.2,
fun.data = median_hilow,
fun.args = list(conf.int = 0.5)
)+
  labs(
    # title = "Model predictive performance", 
    # subtitle = "Model prediction improves with chains but not iterations",
    x = "Number of seeds",
    y = "ARI"
  )


all_results_df %>%
  ggplot(aes(x = N_iter, y = ARI)) +
  geom_line(aes(group = Simulation), colour = "grey", alpha = 0.3) +
  facet_wrap(~N_seeds, labeller = labeller(N_seeds = seed_labels)) + #, ncol = 1) +
  geom_smooth(color = "blue", lty = 1,
              stat = "summary",
              fill = "red",
              alpha = 0.2,
              fun.data = median_hilow,
              fun.args = list(conf.int = 0.5)
  )

p2 <- all_results_df %>%
  ggplot(aes(x = N_seeds, y = Frobenius_product)) +
  geom_line(aes(group = Simulation), colour = "grey", alpha = 0.3) +
  facet_wrap(~N_iter, labeller = labeller(N_iter = iter_labels), ncol = 1) +
  geom_smooth(color = "blue", lty = 1,
              stat = "summary",
              fill = "red",
              alpha = 0.2,
              fun.data = median_hilow,
              fun.args = list(conf.int = 0.5)
  ) +
  labs(
    # title = "Model uncertainty quantification", 
    # subtitle = "Model certainty does not improve with iterations or chains",
    x = "Number of seeds",
    y = "Normalised Frobenius product"
  )

p1+p2 + plot_annotation(
  title = 'Model performance on 2D dataset',
  # subtitle = 'These 3 plots will reveal yet-untold secrets about our beloved data-set',
  subtitle = 'Model predictive performance increases with chains but not iterations; \nUncertainty quantification does not improve with iterations or chains'
)

ggsave("/Users/stephen/Desktop/simple2d_modelperfmance.png")

all_results_df %>%
  ggplot(aes(x = N_iter, y = Frobenius_product)) +
  geom_line(aes(group = Simulation), colour = "grey", alpha = 0.3) +
  facet_wrap(~N_seeds, labeller = labeller(N_seeds = seed_labels)) + #, ncol = 1) +
  geom_smooth(color = "blue", lty = 1,
              stat = "summary",
              fill = "red",
              alpha = 0.2,
              fun.data = median_hilow,
              fun.args = list(conf.int = 0.5)
  )
