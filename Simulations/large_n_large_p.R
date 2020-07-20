
library(mdiHelpR)
library(magrittr)

N <- 3000
P <- 3000

K <- round(log(N))
delta_mu <- c(0.2, 0.4, 0.6)
L <- length(delta_mu)
my_data <- vector("list", L)

for(l in 1:L){
  
  my_data[[l]] <- generateSimulationDataset(K, N, P, delta_mu = delta_mu[l])
}

annotatedHeatmap(my_data[[1]]$data, my_data[[1]]$cluster_IDs)
annotatedHeatmap(my_data[[2]]$data, my_data[[2]]$cluster_IDs)
annotatedHeatmap(my_data[[3]]$data, my_data[[3]]$cluster_IDs, cluster_rows=F, rownames = F)


