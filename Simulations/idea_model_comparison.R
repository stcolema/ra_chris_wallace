#!/usr/env/bin/Rscript

# Idea for analysis of method performance

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
library(mclust)

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


main_dir <- "./PhD/Year_1/Consensus_inference/Consensus_inference_gen/Simulations/Simulation_results/Simulations/Single_dataset/"
scn <- "base_case"
sims_used <- c(1, 2, 4, 6, 8, 10)
n_sim <- length(sims_used)
# sim_num <- 1

n_chains <- c(1, 10, 30, 100)
n_iter <- c(1, 50, 500, 5000)
n_consensus_models <- nrow(expand.grid(n_chains, n_iter))

burnin <- 5e5
thin_by <- 100
expected_length <- 1e6

model_df <- expand.grid(n_chains, n_iter) %>%
  set_colnames(c("n_chains", "n_iter")) %>%
  rbind(c(1, 1e6)) %>%
  rbind(c(NA, NA)) %>%
  cbind(data.frame("Inference" = c(
    rep("Consensus", n_consensus_models),
    "Bayesian",
    "Frequentist"
  )))

model_df$model <- model_df %>%
  with(paste(Inference, n_chains, n_iter, sep = "_"))

model_df$model[17:18] <-  c("Bayesian", "Frequentist")


n_models <- nrow(model_df)

results <- matrix(0, nrow = n_models, ncol = n_sim) %>%
  set_colnames(paste0("Simulation_", sims_used)) %>%
  as.data.frame()

uncertainty_df <- matrix(0, nrow = n_models, ncol = n_sim) %>%
  set_colnames(paste0("Simulation_", sims_used)) %>%
  as.data.frame()


# result_df <- cbind(model_df, results)

# models <- c("Bayesian", paste0("Consensus", ))

for (l in 1:n_sim) {
  sim_num <- sims_used[l]
  truth <- readRDS(paste0(data_dir, "cluster_IDs_", sim_num, ".rds"))
  N <- length(truth)
  truth_coclustering_mat <- cltoSim(truth)
  
  my_data <- read.csv(paste0(data_dir, "dataset_", sim_num, ".csv"), row.names = 1)

  curr_dir <- paste0(main_dir, scn, "/simulation_", sim_num)

  consensus_dir <- paste0(curr_dir, "/Consensus/")
  bayes_dir <- paste0(curr_dir, "/Bayesian/")

  # Read the files and sort by chain number (recogininsing as a numeric)
  consensus_files <- list.files(consensus_dir, full.names = T) %>%
    stringr::str_sort(numeric = T)

  n_consensus <- length(consensus_files)

  consensus_samples <- list()
  for (i in 1:n_consensus) {
    consensus_samples[[i]] <- fread(consensus_files[i], drop = 1)
  }

  # chain_length <- 5
  for (chain_length in n_iter) {
    for (chains_used in n_chains) {
      curr_model_ind <- which(model_df$n_chains == chains_used & model_df$n_iter == chain_length)
      consensus_output <- lapply(consensus_samples[1:chains_used], function(x) {
        x[chain_length, ]
      }) %>%
        unlist() %>%
        matrix(nrow = chains_used, byrow = T)

      consensus_mat <- createSimilarityMat(t(consensus_output))

      consensus_clustering <- mcclust::maxpear(consensus_mat)

      # annotatedHeatmap(scale(my_data), consensus_clustering$cl)
      results[curr_model_ind, l] <- arandi(truth, consensus_clustering$cl)
      
      uncertainty_df[curr_model_ind, l] <- (truth_coclustering_mat - consensus_mat)**2 %>% 
        sum() /  (sum(truth_coclustering_mat - N)) # (N*N)
    }
  }
  # Read the files and sort by chain number (recogininsing as a numeric)
  bayes_files <- list.files(bayes_dir, full.names = T) %>%
    stringr::str_sort(numeric = T)
  n_bayes <- length(bayes_files)

  bayes_samples <- list()
  best_bayes <- 0
  for (i in 1:n_bayes) {
    bayes_samples <- fread(bayes_files[i], drop = 1)

    bayes_psm <- bayes_samples[seq(burnin, expected_length, by = thin_by), ] %>%
      as.matrix() %>%
      t() %>%
      createSimilarityMat()

    # bayes_psm_1 %>%
    #   pheatmap(color = sim_col_pal, breaks = sim_breaks)

    bayes_clustering <- mcclust::maxpear(bayes_psm)

    # annotatedHeatmap(scale(my_data), bayes_clustering$cl)
    bayes_score <- arandi(truth, bayes_clustering$cl)
    if (bayes_score > best_bayes) {
      results[17, l] <- best_bayes <- bayes_score
      
      uncertainty_df[17, l] <- (truth_coclustering_mat - bayes_psm)**2 %>% 
        sum() / (sum(truth_coclustering_mat - N)) # (N*N)
    }
  }

  freq_mod <- mclust::Mclust(my_data, G = 2:30)

  results[18, l] <- arandi(truth, freq_mod$classification)
}

results_df <- cbind(model_df, results) %>%
  tidyr::pivot_longer(dplyr::contains("Simulation"), 
                      names_to = "Simulation", 
                      values_to = "ARI",
                      names_ptypes = list(Simulation = factor())
                      )

uncertainty_df_plt <- cbind(model_df, uncertainty_df) %>%
  tidyr::pivot_longer(dplyr::contains("Simulation"), 
                      names_to = "Simulation", 
                      values_to = "Uncertainty",
                      names_ptypes = list(Simulation = factor())
                      )

# results_df$model <- results_df %>%
#   with(paste(Inference, n_chains, n_iter, sep = "_"))
# 
# results_df$model[17:18] <-  c("Bayesian", "Frequentist")
# results_df$model <- as.factor(results_df$model)
# 
# results_df <- results_df %>%
#   tidyr::pivot_longer(dplyr::contains("Simulation"), names_to = "Simulation", values_to = "ARI", names_ptypes = list(Simulation = factor()))

ggplot(results_df, aes(x = model, y = ARI, group = model)) +
  geom_boxplot(fill = "grey") +
  geom_point(alpha = 0.3) +
  coord_flip() +
  labs(
    title = "Base case",
    subtitle = "Model performance across simulations",
    x = "Model"
  )

ggsave("./PhD/Year_1/Consensus_inference/Consensus_inference_gen/Simulations/Simulation_pipe/example_model_comparison.png")


ggplot(uncertainty_df_plt %>% dplyr::filter(model!="Frequentist"), aes(x = model, y = Uncertainty, group = model)) +
  geom_boxplot(fill = "grey") +
  geom_point(alpha = 0.3) +
  coord_flip() +
  labs(
    title = "Base case",
    subtitle = "Model uncertainty across simulations",
    x = "Model"
  )

ggsave("./PhD/Year_1/Consensus_inference/Consensus_inference_gen/Simulations/Simulation_pipe/example_model_uncertainty_1.png")

# library(mcclust.ext)
# 
# cb_1 <- credibleball(bayes_clustering$cl,bayes_samples)
# 
# 
# 
# truth
# 
# # True allocation for test data
# reference <- factor(truth, levels = unique(truth))
# 
# # Confusion matrix for current fold
# cmlist[[i]] <- conf <- caret::confusionMatrix(
#   data = mcmc_predictions,
#   reference = reference
# )$table
# 
# # Create allocation matrices for truth, filled initially with 0's
# allocmatrix <- matrix(0,
#                       nrow = 200,
#                       ncol = 50
# )
# 
# test_alloc_2 <- allocmatrix
# 
# # The numbers associated with the classes (i.e. numerical representation of
# # the classes)
# class_numerics <- seq(1, length(unique(truth)))
# 
# # create allocation matrix
# for (j in seq_along(1:200)) {
#   # The class the current individual belongs to
#   alloc <- as.numeric(truth, class_numerics)[j]
#   
#   # Enter this in the allocation matrix
#   allocmatrix[j, alloc] <- 1
#   # test_alloc_2[j, map_predictions[j]] <- 1
# }
# 
# bayes_samples[1:100,] %>% t() %>%  unique() %>% length()
# 
# apply(bayes_samples[1:100,],1,unique)
# 
# psm_t <- cltoSim(truth) 
# 
# (bayes_psm - psm_t) %>% sum()
# (consensus_mat - psm_t) %>% sum() / 200**2
#   
#   pheatmap(cluster_rows=F,cluster_cols=F)
#   

rev_truth <- truth %>% rev() %>% cltoSim()

truth_coclustering_mat %>% pheatmap(cluster_rows = F,cluster_cols = F)
rev_truth %>% pheatmap(cluster_rows = F,cluster_cols = F)

comp <- diag(1,nrow=N)
pheatmap(comp)

(truth_coclustering_mat - comp)**2 %>% pheatmap()

(truth_coclustering_mat - comp)**2 %>% sum() / (sum(truth_coclustering_mat)-N)

uncertainty_df_2 <- uncertainty_df * (N**2)/ (sum(truth_coclustering_mat)-N)

uncertainty_df_2_plt <- cbind(model_df, uncertainty_df_2) %>%
  tidyr::pivot_longer(dplyr::contains("Simulation"), 
                      names_to = "Simulation", 
                      values_to = "Uncertainty",
                      names_ptypes = list(Simulation = factor())
  )

ggplot(uncertainty_df_2_plt %>% dplyr::filter(model!="Frequentist"), aes(x = model, y = Uncertainty, group = model)) +
  geom_boxplot(fill = "grey") +
  geom_point(alpha = 0.3) +
  coord_flip() +
  labs(
    title = "Base case",
    subtitle = "Model uncertainty across simulations",
    x = "Model"
  )


ggsave("./PhD/Year_1/Consensus_inference/Consensus_inference_gen/Simulations/Simulation_pipe/example_model_uncertainty_2.png")
