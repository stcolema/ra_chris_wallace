
library(mdiHelpR)
library(magrittr)
library(stringr)

N <- 3000
P <- 3000

K <- round(log(N))
delta_mu <- c(0.2, 0.4, 0.6, 0.8, 1.0)
L <- length(delta_mu)
my_data <- vector("list", L)

my_dir <- "C:/Users/stephen/Documents/PhD/Year_1/Consensus_inference/ra_chris_wallace/Simulations/"

filenames <- paste0(my_dir, "SimulationLargeNLargePDM") %>%
  paste0(str_remove(delta_mu, "\\.")) %>%
  paste0(".png")

title <- "Large N Large P"

for (l in 1:L) {
  main <- paste(title, "DM", delta_mu[l])
  my_data[[l]] <- generateSimulationDataset(K, N, P, delta_mu = delta_mu[l])
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

# annotatedHeatmap(my_data[[1]]$data, my_data[[1]]$cluster_IDs)
# annotatedHeatmap(my_data[[2]]$data, my_data[[2]]$cluster_IDs)
# annotatedHeatmap(my_data[[3]]$data, my_data[[3]]$cluster_IDs, cluster_rows = F, rownames = F)
