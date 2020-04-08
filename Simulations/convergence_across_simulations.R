
library(mdiHelpR)
library(magrittr)
library(coda)
library(mcclust)
library(ggplot2)
library(dplyr)

setMyTheme()

file_dir <- "../Consensus_inference_gen/Analysis/MDI_test_data/Easier_convergence/MATLAB_MDI/N_clust_50/Abbreviated_data/"

seeds <- list.dirs(file_dir, recursive = F)

data_dir <- seeds[stringr::str_order(seeds, decreasing = F, numeric = T)]

list.files(data_dir[1])

N_sim <- 10
N <- 200
N_iter <- 1e4
K <- 8
item_names <- paste0("Person_", 1:N)
param <- "Mass_Parameter_1"
n_model <- 10
clusterings <- list()
models <- list()
g_plots <- list()

truth <- sample(1:K, N, replace = T)

for(l in 1:N_sim){
  models[[l]] <- list()
  clusterings[[l]] <- list()
for(i in 1:n_model){
clusterings[[l]][[i]] <- sample(1:K, size = N*N_iter, replace = T) %>% 
  matrix(nrow = N_iter) %>% 
  set_colnames(item_names)

shape <- rnorm(1, mean = 3)**2
models[[l]][[i]] <- data.frame(MassParameter_1 = rgamma(N_iter, shape, 0.2)) %>% #,
                          # MassParameter_2 = rgamma(N_iter, shape, 0.5)) %>%
  # cbind(clusterings_1) %>%
  mcmc()
}
  g_plots[[l]] <- gelmanPlot(models[[l]])
  rel_data <- g_plots[[l]]$data %>% 
    tibble::add_column(Simulation = l)
  
  if(l > 1){
    gelman_df <- rbind(gelman_df, rel_data)
  } else {
    gelman_df <- rel_data
  }
}


sum_data <- gelman_df %>% 
  group_by(Last_iter) %>% 
  filter(Quantity == "median") %>% 
  summarise(Shrinkage_factor = median(Shrinkage_factor)) %>% 
  tibble::add_column(Quantity = "Summary", Simulation = NA)

p <- gelman_df %>% 
  filter(Quantity == "median") %>% 
  # bind_rows(sum_data)
  # tibble::add_column("Group" = interaction(.$Quantity, .$Simulation)) %>% 
  ggplot(aes(x = Last_iter, y = Shrinkage_factor)) +
  geom_line(aes(group = Simulation, colour = "Simulation"), lty = 3) +
  geom_smooth(aes( color = 'Median'), 
              stat = 'summary',
              fill = 'red', 
              alpha = 0.2, 
              fun.data = median_hilow, 
              fun.args = list(conf.int = 0.5)
            ) +
  # stat_summary(aes(colour = "Median"), fun = "median", geom = "line") +
  labs(
    title = "Convergence across chains",
    subtitle = "Including median and interquartile range",
    x = "Last iteration in chain",
    y = "Gelman-Rubin shrinkage factor",
    colour = "Values"
  ) + 
  ggplot2::geom_hline(aes(yintercept = 1L, colour="Target"), linetype = 2) +
  scale_color_manual(values = c("blue", "grey", "red")) 
  # scale_color_discrete(name = "Values", labels = c("Y2", "Y1"))
  
p

my_tib <- tibble(
  plot_data = list(p$data),
  call = 'gelman_df %>% 
  filter(Quantity == "median") %>% 
  ggplot(aes(x = Last_iter, y = Shrinkage_factor)) +
  geom_line(aes(group = Simulation, colour = "Simulation"), lty = 3) +
  geom_smooth(aes( color = "Median"), 
              stat = "summary",
              fill = "red", 
              alpha = 0.2, 
              fun.data = median_hilow, 
              fun.args = list(conf.int = 0.5)
            ) +
  labs(
    title = "Convergence across chains",
    subtitle = "Including median and interquartile range",
    x = "Last iteration in chain",
    y = "Gelman-Rubin shrinkage factor",
    colour = "Values"
  ) + 
  ggplot2::geom_hline(aes(yintercept = 1L, colour="Target"), linetype = 2) +
  scale_color_manual(values = c("blue", "grey", "red"))'
)

saveRDS(my_tib, file = "~/PhD/Year_1/Consensus_inference/Consensus_inference_gen/Simulations/Simulation_pipe/convergence_plot_tib.Rds")

# tib1 <- readRDS("~/PhD/Year_1/Consensus_inference/Consensus_inference_gen/Simulations/Simulation_pipe/convergence_plot_tib.Rds")

ggsave("~/PhD/Year_1/Consensus_inference/Consensus_inference_gen/Simulations/Simulation_pipe/example_converge_across_sims_plot_band.png")


sim_1 <- models[[1]]

sim_1[[1]] %>% 
  as.data.frame() %>% 
  add_column(Simulation = "Simulation 1") %>% 
  head()

apply(clusterings[[1]][[1]], 1, arandi, truth)

psm1 <- comp.psm(clusterings[[1]][[1]])

pred_clustering <- maxpear(psm1)

cb <- mcclust.ext::credibleball(pred_clustering$cl, clusterings[[1]][[1]])

cb$c.uppervert %>% 
  arandi(pred_clustering$cl)

cb$c.uppervert %>% 
  arandi(truth)


cb$c.lowervert %>% 
  arandi(pred_clustering$cl)

cb$c.lowervert %>% c() %>% 
  unique()

  arandi(truth)

cl_df <- data.frame(
  "Predicted" = pred_clustering$cl,
  "Upper" = c(cb$c.uppervert),
  "Lower" = c(cb$c.lowervert)
  )


for(i in 1:10){
  sim_df <- clusterings[[1]][[i]] %>% 
    as.data.frame() %>% 
    add_column(Simulation = paste0("Simulation ", i)) 
  if(i > 1){
    uncertainty_df <- bind_rows(uncertainty_df, sim_df)
  } else {
    uncertainty_df <- sim_df
  }
}

uncertainty_df$Simulation <- uncertainty_df$Simulation %>% 
  factor(levels = stringr::str_sort(unique(.), numeric = T))

# ggplot(uncertainty_df, )





