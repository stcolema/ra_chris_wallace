#!/usr/env/bin/Rscript

# Various libraries
library(pheatmap)
library(ggplot2)
library(mdiHelpR)
library(magrittr)
library(coda)
library(mcclust)
library(data.table)
library(patchwork)
library(stringr)

# Not used
# library(fitR) # devtools::install_github("sbfnk/fitR")

# My dewfault ggplot2 theme
setMyTheme()

# Colour and breaks for PSMs
sim_col_pal <- mdiHelpR::simColPal()
sim_breaks <- defineBreaks(sim_col_pal, lb = 0)

col_pal <- dataColPal()
cor_breaks <- defineBreaks(col_pal)

# Directory
data_dir <- "./PhD/Year_1/Consensus_inference/Consensus_inference_gen/Simulations/Input data/Datasets/base_case/"
my_data <- read.csv(paste0(data_dir, "dataset_1.csv"), row.names = 1)
truth <- readRDS(paste0(data_dir, "cluster_IDs_1.rds"))

main_dir <- "./PhD/Year_1/Consensus_inference/Consensus_inference_gen/Simulations/Simulation_results/Simulations/Single_dataset/"
scn <- "base_case"
sim_num <- 1
curr_dir <- paste0(main_dir, scn, "/simulation_", sim_num)

consensus_dir <- paste0(curr_dir, "/Consensus/")
bayes_dir     <- paste0(curr_dir, "/Bayesian/")

# Read the files and sort by chain number (recogininsing as a numeric)
conesnsus_files <- list.files(consensus_dir, full.names = T) %>%
  stringr::str_sort(numeric = T)

n_consensus <- length(conesnsus_files)
chain_length <- 5
consensus_samples <- list()
for (i in 1:n_consensus) {
  consensus_samples[[i]] <- fread(conesnsus_files[i], drop = 1)
}

consensus_output <- lapply(consensus_samples, function(x){ x[chain_length, ]}) %>% 
  unlist() %>% 
  matrix(nrow = 100, byrow = T)

# Quick check that the different seeds are behaving correctly
pheatmap(consensus_output, cluster_rows = F)


consensus_output %>% 
  cor() %>% 
  pheatmap(color = col_pal, breaks = cor_breaks)

# Check correlation of chains
consensus_output %>% 
  t() %>% 
  cor() %>% 
  pheatmap(color = col_pal, breaks = cor_breaks)

consensus_psm <- createSimilarityMat(t(consensus_output))

consensus_psm %>% 
  pheatmap(color = sim_col_pal, breaks = sim_breaks)

consensus_clustering <- mcclust::maxpear(consensus_psm)

annotatedHeatmap(scale(my_data), consensus_clustering$cl)
arandi(truth, consensus_clustering$cl)

# Read the files and sort by chain number (recogininsing as a numeric)
bayes_files <- list.files(bayes_dir, full.names = T) %>%
  stringr::str_sort(numeric = T)
n_bayes <- length(bayes_files)

bayes_samples <- list()
for (i in 1:n_bayes) {
  bayes_samples[[i]] <- fread(bayes_files[i], drop = 1)
}

bayes_psm_1 <- bayes_samples[[1]][seq(5e5, 6e5, by = 100), ] %>% 
  as.matrix() %>% 
  t() %>% 
  createSimilarityMat()

bayes_psm_1 %>% 
  pheatmap(color = sim_col_pal, breaks = sim_breaks)

bayes_clustering <- mcclust::maxpear(bayes_psm_1)

annotatedHeatmap(scale(my_data), bayes_clustering$cl)
arandi(truth, bayes_clustering$cl)

# 
# 
# gaussianBIC <- function(x, k, ll, diagonal_cov = F){
#   
#   p <- ncol(x)
#   n <- nrow(x)
#   
#   # For each mixture:
#   # Contribution of covariance matrix to parameter count
#   cov_params <- ((p*p - p)/2 + p) * (1 - diagonal_cov) + p * diagonal_cov
#   
#   # Contribution of mean vector
#   mean_params <- p
#   
#   # Mixture weight
#   wgt_param <- 1
#   
#   # Symmetrix covariance matrix (p*p - p) + p, mean vector p and mixture weight 1
#   mixture_params <- cov_params + mean_params + wgt_param
#   
#   # We have K mixtures but the weights are constrained to sum to 1
#   n_param <- k * mixture_params - 1
#   
#   n_param * log(n) - 2 * ll
# }