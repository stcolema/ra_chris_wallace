#!/usr/bin/env Rscript

# === Libraries ================================================================

library(magrittr)
library(rjags)
library(fitR) # some extensions for rjags
library(ggplot2)

# === Yeast data with some generated clusters ==================================

yeast_data_dir <- "~/Desktop/Yeast_3_datasets/MDI_output"

# Directroy to save plots to
yeast_plot_dir <- "~/Desktop/ra_chris_wallace/Notes/Thesis/Images/Gen_data/Case_1"


all_dirs <- list.dirs(yeast_data_dir, recursive = F)
all_types <- list.dirs(yeast_data_dir, recursive = F, full.names = F)

burn <- 20000 

long_dirs <- all_dirs %>% 
  extract(grep("Long", .))

long_types <- all_types %>% 
  extract(grep("Long", .))

long_seed <- long_types %>% 
  regmatches(., gregexpr("[[:digit:]]+", .)) %>%
  unlist()

consensus_dirs <- all_dirs %>% 
  extract(grep("Consensus", .))

consensus_types <- all_types %>% 
  extract(grep("Consensus", .))

consensus_iter <- consensus_types %>% 
  regmatches(., gregexpr("[[:digit:]]+", .)) %>%
  unlist()

consensus_files <- paste0("consensus_",  consensus_iter, ".csv") %>% 
  paste0(consensus_dirs, "/", .)

long_files <- paste0("out_seed_" , long_seed, ".csv") %>% 
  paste0(long_dirs, "/", .)

n_long <- length(long_files)
n_consensus <- length(consensus_files)
n_param <- 6

long_mcmc_yeast <- list()
consensus_mcmc_yeast <- list()

for(i in 1:n_long){
  curr_data <- read.csv(long_files[i], header = T)[1:40000, 1:n_param] # stop at 40,000 as 2 million is good (and differing lengths)
  
  long_mcmc_yeast[[i]] <- mcmc(curr_data)
}

# 
for(i in 1:n_consensus){
  curr_data <- read.csv(consensus_files[i], header = T)[, 1:6]
  
  consensus_mcmc_yeast[[i]] <- mcmc(curr_data)
  
}




# gelman.plot(consensus_mcmc_yeast)
png(paste0(yeast_plot_dir, "/Gelman_plot.png"))
gelman.plot(long_mcmc_yeast)
dev.off()

for(i in 1:n_long){
  png(paste0(yeast_plot_dir, "/Esimated_burn_in_plot_", long_seed[i], ".png"))
  plotESSBurn(long_mcmc_yeast[[i]])
  dev.off()
  
  png(paste0(yeast_plot_dir, "/Auto_correlation_plot_", long_seed[i], ".png"))
  autocorr.plot(long_mcmc_yeast[[i]])
  dev.off()
}

burn_mcmc_yeast <- list()

# Burn in for data
for(i in 1:n_long){
  # curr_data <- read.csv(long_files[i], header = T)[10000:40000, 1:6] # stop at 40,000 as 2 million is good (and differing lengths)
  
  burn_mcmc_yeast[[i]] <- burnAndThin(long_mcmc_yeast[[i]], burn = 10000, thin = 1)
}

# Individual chain stationary
for(i in 1:n_long){
  png(paste0(yeast_plot_dir, "/Geweke_plot_burn_", burn, "_", long_seed[i], ".png"))
  geweke.plot(burn_mcmc_yeast[[i]])
  dev.off()
}


p_values_yeast <- vector("list", length = n_long)
z_scores_yeast <- vector("list", length = n_long)
for(i in 1:n_long){
  gew_diag <- geweke.diag(burn_mcmc_yeast[[i]])
  
  p_values_yeast[[i]] <- rep(0, n_param)
  names(p_values_yeast[[i]]) <- names(gew_diag$z)
  
  z_scores_yeast[[i]] <- gew_diag$z
  
  for(j in 1:n_param){
    # Ca;culate a p-value from the z-score
    p_value <- 2*pnorm(-abs(gew_diag$z[j]))
    p_values_yeast[[i]][j] <- p_value
  }
}

p_values_df_yeast <- data.frame(do.call(rbind, p_values_yeast))
z_values_df_yeast <- data.frame(do.call(rbind, z_scores_yeast))

# === Generated data ===========================================================

gen_data_dir <- "~/Desktop/Gen_data_output/MDI_output"


# Directroy to save plots to
gen_plot_dir <- "~/Desktop/ra_chris_wallace/Notes/Thesis/Images/Gen_data/Case_2"
dir.create(plot_dir, showWarnings = F)

all_dirs <- list.dirs(gen_data_dir, recursive = F)
all_types <- list.dirs(gen_data_dir, recursive = F, full.names = F)

long_dirs <- all_dirs %>% 
  extract(grep("Long", .))

long_types <- all_types %>% 
  extract(grep("Long", .))

long_seed <- long_types %>% 
  regmatches(., gregexpr("[[:digit:]]+", .)) %>%
  unlist()

consensus_dirs <- all_dirs %>% 
  extract(grep("Consensus", .))

consensus_types <- all_types %>% 
  extract(grep("Consensus", .))

consensus_iter <- consensus_types %>% 
  regmatches(., gregexpr("[[:digit:]]+", .)) %>%
  unlist()

consensus_files <- paste0("gen_data_",  consensus_iter, "_iter_1000.csv") %>% 
  paste0(consensus_dirs, "/", .)

long_files <- paste0("out_seed_" , long_seed, ".csv") %>% 
  paste0(long_dirs, "/", .)

n_long <- length(long_files)
n_consensus <- length(consensus_files)
n_param <- 6

long_mcmc_gen <- list()
consensus_mcmc_gen <- list()

for(i in 1:n_long){
    curr_data <- read.csv(long_files[i], header = T)[1:40000, 1:6] # stop at 40,000 as 2 million is good (and differing lengths)
    
    long_mcmc_gen[[i]] <- mcmc(curr_data)
}

# 
for(i in 1:n_consensus){
  curr_data <- read.csv(consensus_files[i], header = T)[, 1:n_param]

  consensus_mcmc_gen[[i]] <- mcmc(curr_data)
}
# gelman.plot(consensus_mcmc_gen)

png(paste0(gen_plot_dir, "/Gelman_plot.png"))
gelman.plot(long_mcmc_gen)
dev.off()


# ggsave(paste0(gen_plot_dir, "/Gelman_plot.png"))


for(i in 1:n_long){
  png(paste0(gen_plot_dir, "/Esimated_burn_in_plot_", long_seed[i], ".png"))
  plotESSBurn(long_mcmc_gen[[i]])
  dev.off()
  
  png(paste0(gen_plot_dir, "/Auto_correlation_plot_", long_seed[i], ".png"))
  autocorr.plot(long_mcmc_gen[[i]])
  dev.off()
}

burn_mcmc_gen <- list()

# Burn in for data
for(i in 1:n_long){
  # curr_data <- read.csv(long_files[i], header = T)[10000:40000, 1:6] # stop at 40,000 as 2 million is good (and differing lengths)
  
  burn_mcmc_gen[[i]] <- burnAndThin(long_mcmc_gen[[i]], burn = burn, thin = 1)
}

# Individual chain stationary
for(i in 1:n_long){
  png(paste0(gen_plot_dir, "/Geweke_plot_burn_", burn, "_", long_seed[i], ".png"))
  geweke.plot(burn_mcmc_gen[[i]])
  dev.off()
  
  png(paste0(gen_plot_dir, "/Auto_correlation_plot_burn_", burn, "_", long_seed[i], ".png"))
  autocorr.plot(burn_mcmc_gen[[i]])
  dev.off()
}

p_values_gen <- vector("list", length = n_long)
z_scores_gen <- vector("list", length = n_long)
for(i in 1:n_long){
  gew_diag <- geweke.diag(burn_mcmc_gen[[i]])
  
  p_values_gen[[i]] <- rep(0, n_param)
  names(p_values_gen[[i]]) <- names(gew_diag$z)
  
  z_scores_gen[[i]] <- gew_diag$z
  
  for(j in 1:n_param){
    # Ca;culate a p-value from the z-score
    p_value <- 2*pnorm(-abs(gew_diag$z[j]))
    p_values_gen[[i]][j] <- p_value
  }
}

p_values_df_gen <- data.frame(do.call(rbind, p_values_gen))
z_values_df_gen <- data.frame(do.call(rbind, z_scores_gen))
gelman.diag(burn_mcmc_gen)

png(paste0(gen_plot_dir, "/Gelman_plot_burn_", burn, ".png"))
gelman.plot(burn_mcmc_gen)
dev.off()

# == P-values ==================================================================

p_values_df_yeast$Data <- "Case 1"
p_values_df_gen$Data <- "Case 2"

p_values_df <- rbind(p_values_df_yeast, p_values_df_gen)

p_values_vec <- c(unlist(p_values_yeast), unlist(p_values_gen))

data_type_vec <- c(rep("Case 1", nrow(p_values_df_yeast)),
                       rep("Case 2", nrow(p_values_df_gen)  )
)

# p_adj_vec <- p.adjust(p_values_vec, method = "BH")

# Put them all in a data.frame
p_adj_mat <- t(matrix(p_values_vec, nrow = 6)) %>% 
  set_colnames(names(gew_diag$z)) %>% 
  round(., 3) %>% 
  as.data.frame()

p_adj_mat$Case <- data_type_vec

# === Z-scores ================================================================

z_values_df <- rbind(z_values_df_yeast, z_values_df_gen)

z_values_vec <- c(unlist(z_scores_yeast), unlist(z_scores_gen))

data_type_vec <- c(rep("Case 1", nrow(z_values_df_yeast)),
                   rep("Case 2", nrow(z_values_df_gen)  )
)

# z_adj_vec <- p.adjust(z_values_vec, method = "BH")

# Put them all in a data.frame
z_mat <- t(matrix(z_values_vec, nrow = 6)) %>% 
  set_colnames(names(gew_diag$z)) %>% 
  round(., 3)

z_df <- z_mat %>% 
  as.data.frame()

z_values_df_yeast[rowSums(abs(z_values_df_yeast) > 1.96) > 0,]
z_values_df_gen[rowSums(abs(z_values_df_gen) > 1.96) > 0,]
# (abs(z_values_df_gen) > 1.96)

z_mat[abs(z_mat) > 1.96]


z_dt$Case <- data_type_vec

# === Nonsense =================================================================

# N <- 1000
# x <- 1:N
# epsilon <- rnorm(N, 0, 1)
# y <- x + epsilon
# 
# jags <- jags.model('/home/MINTS/sdc56/Desktop/rjags_stuff/example.bug',
#                    data = list('x' = x,
#                                'y' = y,
#                                'N' = N),
#                    n.chains = 4,
#                    n.adapt = 10)
# samples <- coda.samples(jags,
#                         c('a', 'b'),
#                         1000)
# png("plot_1.png")
# plot(samples)
# dev.off()
# 
# data(mcmc)
# 
# mcmc_data <- list()
# main_path <- "/home/MINTS/sdc56/Desktop/Yeast_MDI/Yeast_output/Diffuse"
# 
# 
# 
# weight <- "weight_2"
# main_dirs <- list.dirs(path = main_path, full.names = TRUE, recursive = F) %>% 
#   extract(grepl(weight, .))
# 
# long_sub_path <- "/Long_runs"
# long_dirs <- list.dirs(path = paste0(main_path, long_sub_path), full.names = TRUE, recursive = F) %>% 
#   extract(grepl(weight, .))
# 
# all_dirs <- c(main_dirs, long_dirs)
# 
# mcmc_files <- list.files(all_dirs, pattern = "*.csv", full.names = T)
# 
# consensus_files <- mcmc_files[1:2]
# long_files <- mcmc_files[-c(1:2)]
# 
# mcmc_data <- list()
# mcmc_data_short <- list()
# n_files <- length(mcmc_files)
# for(i in 1:n_files){
#   curr_data <- read.csv(mcmc_files[i])[,1:21]
#   
#   mcmc_data[[i]] <- mcmc(curr_data)
#   
#   # mcmc_data_short[[i]] <- mcmc(curr_data[1:100,])
# }

# second_path <- "/home/MINTS/sdc56/Desktop/New_long_06/Diffuse_long_runs_0.6"
# new_data <- list.files(second_path, full.names = T)
# 
# new_mcmc <- list()
# n_2 <- length(new_data)
# for(i in 1:n_2){
#   curr_data <- read.csv(new_data[i])[20000:200000, 1:6]
#   
#   new_mcmc[[i]] <- mcmc(curr_data)
# }
# 
# gelman.plot(new_mcmc)

# mcmc_data %>% str()
# 
# long_mcmc <- mcmc_data[-c(1,2)]
# gelman.plot(long_mcmc)
# 
# consensus_mcmc <- mcmc_data[1:2]
# 
# 
# gelman.plot(consensus_mcmc)
# 
# mcmc_x <- mcmc(x_1)
# summary(mcmc_x)
# acceptanceRate <- 1 - rejectionRate(mcmc_x)
# acceptanceRate
# 
# plot(mcmc_x)
# mcmc.trace <- mcmc_x
# mcmc.trace.burned <- burnAndThin(mcmc_x, burn = 20, thin = 1)
# plot(mcmc.trace.burned)
# autocorr.plot(mcmc.trace.burned)
# plotESSBurn(mcmc.trace)
# 
# 
# gelman.plot(mcmc.trace)
# 
# # samples <- coda.samples(mcmc_data,
# #                         c('MassParameter_2', 'MassParameter_1'),
# #                         10)   
# 
# summary(mcmc_data[[1]])
# plot(mcmc_data[[1]])
# mcmc.trace.burned <- burnAndThin(mcmc_data[[1]], burn = 2900, thin = 1)
# plot(mcmc.trace.burned)
# autocorr.plot(mcmc.trace.burned)
# plotESSBurn(mcmc_data[[1]])
# plotESSBurn(mcmc.trace.burned)
# 
# mcmc.trace.burned_2 <- burnAndThin(mcmc_data[[2]], burn = 5900, thin = 1)
# plot(mcmc.trace.burned_2)
# autocorr.plot(mcmc.trace.burned_2)
# plotESSBurn(mcmc_data[[2]])
# plotESSBurn(mcmc.trace.burned_2)
# 
# plotESSBurn(mcmc_data[[1]]) # burn in to 2900
# plotESSBurn(mcmc_data[[2]]) # burn in to 5900
# plotESSBurn(mcmc_data[[3]]) # no burn in required
# plotESSBurn(mcmc_data[[4]]) # does this converge?
# plotESSBurn(mcmc_data[[5]]) # as 1
# plotESSBurn(mcmc_data[[6]]) # similar to 1 (slightly different)
# plotESSBurn(mcmc_data[[7]]) # between 1 & 3
# plotESSBurn(mcmc_data[[8]]) # as 3 
# plotESSBurn(mcmc_data[[9]]) # as 4
# plotESSBurn(mcmc_data[[10]]) # as 3
# 
# 
# geweke.plot(mcmc_data[[1]])
# geweke.plot(mcmc_data[[2]])
# geweke.plot(mcmc_data[[3]])
# geweke.plot(mcmc_data[[4]])
# geweke.plot(mcmc_data[[5]])
# geweke.plot(mcmc_data[[6]])
# geweke.plot(mcmc_data[[7]])
# geweke.plot(mcmc_data[[8]])
# geweke.plot(mcmc_data[[9]])
# geweke.plot(mcmc_data[[10]])
# 
# mcmc_burned <- list()
# for(i in 1:10){
#   mcmc_burned[[i]] <- burnAndThin(mcmc_data[[i]], burn = 6000, thin = 1)
#   
# }
# 
# geweke.plot(mcmc_burned[[1]])
# geweke.plot(mcmc_burned[[2]])
# geweke.plot(mcmc_burned[[3]])
# geweke.plot(mcmc_burned[[4]])
# geweke.plot(mcmc_burned[[5]])
# geweke.plot(mcmc_burned[[6]])
# geweke.plot(mcmc_burned[[7]])
# geweke.plot(mcmc_burned[[8]])
# geweke.plot(mcmc_burned[[9]])
# geweke.plot(mcmc_burned[[10]])
# 
# 
# 
# gelman.plot(mcmc_burned)
# # gelman.plot(mcmc_data_short)
# 
# geweke.diag(mcmc_data[[1]], frac1=0.1, frac2=0.5)
# geweke.plot(mcmc_data[[1]], frac1=0.1, frac2=0.5)
# 
# summary(many_seeds_mcmc_data[[1]])
# plot(many_seeds_mcmc_data[[1]])
# 
# autocorr.plot(mcmc.trace.burned)
# plotESSBurn(many_seeds_mcmc_data[[1]])
# mcmc.trace.burned <- burnAndThin(many_seeds_mcmc_data[[1]], burn = 0, thin = 1)
# plot(mcmc.trace.burned)
# gelman.plot(many_seeds_mcmc_data)
# 
# # geweke.diag(many_seeds_mcmc_data[[1]], frac1=0.1, frac2=0.5)
# geweke.plot(many_seeds_mcmc_data[[1]], frac1=0.1, frac2=0.5)
# geweke.plot(mcmc.trace.burned)
