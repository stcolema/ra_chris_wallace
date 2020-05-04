
# devtools::install_github("stcolema/mdiHelpR")
library(mdiHelpR)

library(magrittr)
library(pheatmap)
library(mcclust)
library(mclust)

# Where the data lives
data_f <- "/Users/stephen/Desktop/Testing_pipeline/simple_2d/Input_data/dataset_1.csv"
truth_f <- "/Users/stephen/Desktop/Testing_pipeline/simple_2d/Input_data/cluster_IDs_1.Rds"
samples_f <- "/Users/stephen/Desktop/Testing_pipeline/simple_2d/MDI_output/simulation_1/Bayesian/BayesianN1000001T1000Seed1.csv"
my_data <- read.csv(data_f, row.names = 1)
truth <- readRDS(truth_f)

mdiHelpR::annotatedHeatmap(my_data, truth)

samples <- data.table::fread(samples_f, drop = 1)

pheatmap(samples, cluster_rows = F, cluster_cols = F, main = "Samples recorded")

as.matrix(samples[1:4, 1:3], nrow = 4, byrow = T)

dim(samples)

score <- c()
N <- nrow(samples)
samples_used <- c(1, 10, 100, 1000)
n_samples_used <- length(samples_used)

for(i in 1:n_samples_used){
  n_samples <- samples_used[i]
 start_iter <- N - n_samples + 1 
psm <- createSimilarityMat(as.matrix(samples[1:n_samples, ], ncol = 100, byrow = T))

pred_cl <- maxpear(psm)$cl
curr_score <- arandi(truth, pred_cl)
score <- c(score, curr_score)

my_df <- data.frame(ARI = curr_score,
                    Start_iter = 1,
                    End_iter = n_samples,
                    Total_samples = length(1:n_samples)
)

if(i == 1){
  df_of_i <- my_df
} else{
  df_of_i <- rbind(df_of_i, my_df)
}
  # n_samples_used <- c(n_samples_used, length(1:n_samples))
psm <- createSimilarityMat(as.matrix(samples[start_iter:N, ], ncol = 100, byrow = T))
pred_cl <- maxpear(psm)$cl

curr_score <- arandi(truth, pred_cl)
score <- c(score, curr_score)

my_df <- data.frame(ARI = curr_score,
                    Start_iter = start_iter,
                    End_iter = N,
                    Total_samples = length(start_iter:N)
)
  df_of_i <- rbind(df_of_i, my_df)
}

sim_col <- mdiHelpR::simColPal()
my_breaks <- defineBreaks(sim_col, lb =0 )

pheatmap(psm, color = sim_col, breaks = my_breaks)

pheatmap(psm, main = "PSM")

# Try mclust with a subset used in the initialisation
bad_mclust <- Mclust(my_data, G = c(2,5:20), initialization = list(subset = 1:25))
bad_mclust$G
arandi(bad_mclust$classification, truth)

mdiHelpR::annotatedHeatmap(my_data, pred_cl)
mdiHelpR::annotatedHeatmap(my_data, bad_mclust$classification)
