
library(magrittr)
library(rjags)
library(fitR)
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

mcmc_data <- list()
mcmc_data_short <- list()
for(i in 1:10){
  curr_data <- read.csv(paste0("/home/MINTS/sdc56/Desktop/Gen_data_output/Longer_", 
                                    i,
                                    "/out_seed_",
                                    i,
                                    ".csv")
  )[1:40000,1:6]
  
  mcmc_data[[i]] <- mcmc(curr_data)
  
  mcmc_data_short[[i]] <- mcmc(curr_data[1:100,])
}

many_seeds_mcmc_data <- list()
many_seeds_mcmc_data_short <- list()
seeds <- c(
  500,
  1000,
  5000,
  10000
)
n_seeds <- length(seeds)
n_iter <- 1000
for(i in 1:n_seeds){
  curr_seed <- seeds[i]
  curr_data <- read.csv(paste0("/home/MINTS/sdc56/Desktop/Gen_data_output/Many_seeds_", 
                               curr_seed,
                               "/gen_data_",
                               curr_seed,
                                "_iter_1000.csv"
                               )
  )[,1:6]
  
  many_seeds_mcmc_data[[i]] <- mcmc(curr_data)
  
  many_seeds_mcmc_data_short[[i]] <- mcmc(curr_data[1:100,])
}



# x <- read.csv("/home/MINTS/sdc56/Desktop/Gen_data_output/Long_1/out_seed_1.csv")
# x_1 <- x[,1:6]
# head(x_1)


mcmc_x <- mcmc(x_1)
summary(mcmc_x)
acceptanceRate <- 1 - rejectionRate(mcmc_x)
acceptanceRate

plot(mcmc_x)
mcmc.trace <- mcmc_x
mcmc.trace.burned <- burnAndThin(mcmc_x, burn = 20, thin = 1)
plot(mcmc.trace.burned)
autocorr.plot(mcmc.trace.burned)
plotESSBurn(mcmc.trace)


gelman.plot(mcmc.trace)

# samples <- coda.samples(mcmc_data,
#                         c('MassParameter_2', 'MassParameter_1'),
#                         10)   

summary(mcmc_data[[1]])
plot(mcmc_data[[1]])
mcmc.trace.burned <- burnAndThin(mcmc_data[[1]], burn = 2900, thin = 1)
plot(mcmc.trace.burned)
autocorr.plot(mcmc.trace.burned)
plotESSBurn(mcmc_data[[1]])
plotESSBurn(mcmc.trace.burned)

mcmc.trace.burned_2 <- burnAndThin(mcmc_data[[2]], burn = 5900, thin = 1)
plot(mcmc.trace.burned_2)
autocorr.plot(mcmc.trace.burned_2)
plotESSBurn(mcmc_data[[2]])
plotESSBurn(mcmc.trace.burned_2)

plotESSBurn(mcmc_data[[1]]) # burn in to 2900
plotESSBurn(mcmc_data[[2]]) # burn in to 5900
plotESSBurn(mcmc_data[[3]]) # no burn in required
plotESSBurn(mcmc_data[[4]]) # does this converge?
plotESSBurn(mcmc_data[[5]]) # as 1
plotESSBurn(mcmc_data[[6]]) # similar to 1 (slightly different)
plotESSBurn(mcmc_data[[7]]) # between 1 & 3
plotESSBurn(mcmc_data[[8]]) # as 3 
plotESSBurn(mcmc_data[[9]]) # as 4
plotESSBurn(mcmc_data[[10]]) # as 3


geweke.plot(mcmc_data[[1]])
geweke.plot(mcmc_data[[2]])
geweke.plot(mcmc_data[[3]])
geweke.plot(mcmc_data[[4]])
geweke.plot(mcmc_data[[5]])
geweke.plot(mcmc_data[[6]])
geweke.plot(mcmc_data[[7]])
geweke.plot(mcmc_data[[8]])
geweke.plot(mcmc_data[[9]])
geweke.plot(mcmc_data[[10]])

mcmc_burned <- list()
for(i in 1:10){
  mcmc_burned[[i]] <- burnAndThin(mcmc_data[[i]], burn = 6000, thin = 1)
  
}

geweke.plot(mcmc_burned[[1]])
geweke.plot(mcmc_burned[[2]])
geweke.plot(mcmc_burned[[3]])
geweke.plot(mcmc_burned[[4]])
geweke.plot(mcmc_burned[[5]])
geweke.plot(mcmc_burned[[6]])
geweke.plot(mcmc_burned[[7]])
geweke.plot(mcmc_burned[[8]])
geweke.plot(mcmc_burned[[9]])
geweke.plot(mcmc_burned[[10]])


gelman.plot(mcmc_data)
gelman.plot(mcmc_burned)
# gelman.plot(mcmc_data_short)

geweke.diag(mcmc_data[[1]], frac1=0.1, frac2=0.5)
geweke.plot(mcmc_data[[1]], frac1=0.1, frac2=0.5)

summary(many_seeds_mcmc_data[[1]])
plot(many_seeds_mcmc_data[[1]])

autocorr.plot(mcmc.trace.burned)
plotESSBurn(many_seeds_mcmc_data[[1]])
mcmc.trace.burned <- burnAndThin(many_seeds_mcmc_data[[1]], burn = 0, thin = 1)
plot(mcmc.trace.burned)
gelman.plot(many_seeds_mcmc_data)

# geweke.diag(many_seeds_mcmc_data[[1]], frac1=0.1, frac2=0.5)
geweke.plot(many_seeds_mcmc_data[[1]], frac1=0.1, frac2=0.5)
geweke.plot(mcmc.trace.burned)
