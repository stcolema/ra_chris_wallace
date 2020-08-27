
library(mdiHelpR)
library(magrittr)
library(stringr)

set.seed(1)

N <- 200
P <- 20

K <- round(log(N))
pi <- c(0.5, 0.25, 0.125, 0.0625, 0.0625)
delta_mu <- c(0.2, 0.4, 1)
L <- length(delta_mu)
my_data <- vector("list", L)

my_dir <- "C:/Users/stephen/Documents/PhD/Year_1/Consensus_inference/ra_chris_wallace/Simulations/"

filenames <- paste0(my_dir, "SimulationVaryingProportionsDM") %>%
  paste0(str_remove(delta_mu, "\\.")) %>%
  paste0(".png")

titles <- paste("Varying proportions, dm", delta_mu)

for (l in 1:L) {
  main <- titles[l]
  my_data[[l]] <- generateSimulationDataset(K, N, P, delta_mu = delta_mu[l], pi = pi)
  col_order <- mdiHelpR::findOrder(t(my_data[[l]]$data))
  
  annotatedHeatmap(my_data[[l]]$data[, col_order],
                   my_data[[l]]$cluster_IDs,
                   cluster_cols = F,
                   show_rownames = F,
                   show_colnames = F,
                   main = main,
                   silent = T,
                   filename = filenames[l]
  )
}
