

library(mdiHelpR)
library(magrittr)
set.seed(1)

P <- N <- 3e3
K <- round(log(N))
dm <- c(1, 0.5, 0.3)

my_data_1 <- generateSimulationDataset(K, N, P, delta_mu = dm[1])
my_data_2 <- generateSimulationDataset(K, N, P, delta_mu = dm[2])
my_data_3 <- generateSimulationDataset(K, N, P, delta_mu = dm[3])

annotatedHeatmap(my_data_1$data, my_data_1$cluster_IDs, 
                 show_rownames = F,
                 show_colnames = F, 
                 main = "N = 3k, P = 3k, dm = 1",
                 filename= "../LargeNLargePDm10.png")

annotatedHeatmap(my_data_2$data, my_data_2$cluster_IDs, 
                 show_rownames = F,
                 show_colnames = F,
                 main = "N = 3k, P = 3k, dm = 0.5",
                 filename= "../LargeNLargePDm05.png")


annotatedHeatmap(my_data_3$data, my_data_3$cluster_IDs, 
                 show_rownames = F,
                 show_colnames = F,
                 main = "N = 3k, P = 3k, dm = 0.3",
                 filename= "../LargeNLargePDm03.png")
