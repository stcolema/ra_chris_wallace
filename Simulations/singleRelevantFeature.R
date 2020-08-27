
library(mdiHelpR)
library(magrittr)
library(stringr)

set.seed(1)

N <- 200
P <- 1
P_n <- c(1, 100, 1000, 3000)

K <- round(log(N))
delta_mu <- 1
L <- length(P_n)
my_data <- vector("list", L)

my_dir <- "C:/Users/stephen/Documents/PhD/Year_1/Consensus_inference/ra_chris_wallace/Simulations/"

filenames <- paste0(my_dir, "SimulationSingleFeaturePn") %>%
  paste0(P_n) %>%
  paste0(".png")

titles <- paste("One relevant feature,", P_n, "irrelevant features")

for (l in 1:L) {
  main <- titles[l]
  my_data[[l]] <- generateSimulationDataset(K, N, P, p_n = P_n[l])
  if(P_n[l] > 1){
  col_order <- mdiHelpR::findOrder(t(my_data[[l]]$data))
  my_data[[l]]$data <- my_data[[l]]$data[, col_order]
  }
  
  annotatedHeatmap(my_data[[l]]$data,
                   my_data[[l]]$cluster_IDs,
                   cluster_cols = F,
                   show_rownames = F,
                   show_colnames = F,
                   main = main,
                   silent = T,
                   filename = filenames[l]
  )
}
