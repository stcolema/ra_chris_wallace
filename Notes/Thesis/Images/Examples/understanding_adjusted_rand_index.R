#!/usr/bin/env Rscript

set.seed(1)
K <- 3
n <- 1e4

n_1 <- (7/16) * n
n_2 <- (9/16) * n

x <- c(rep(1, n_1), sample(1:K, size = n_2, replace = T))
y <- c(rep(1, n_1), sample(1:K, size = n_2, replace = T))
(cont_table_2 <- table(x, y) / n ) 

mcclust::arandi(x , y)
fossil::rand.index(x, y)

x_2 <- sample(1:K, size = n, replace = T)
y_2 <- sample(1:K, size = n, replace = T)
(cont_table_2 <- table(x_2, y_2) / n ) 

mcclust::arandi(x_2 , y_2)
fossil::rand.index(x_2, y_2)

# Data_x_y <- data.frame(X = x_2, Y = y_2)
# 
# for(i in 1:K){
# Data_x_y$X[Data_x_y$X == i & Data_x_y$Y == i] <- (i %% K) + 1
# Data_x_y$X[Data_x_y$X == i & Data_x_y$Y == i] <- (i %% K) + 2
# }
# 
# (cont_table_2 <- table(Data_x_y$X, Data_x_y$Y) / n ) 
# 
# 
# mcclust::arandi(Data_x_y$X , Data_x_y$X)
# fossil::rand.index(Data_x_y$X, Data_x_y$X)

