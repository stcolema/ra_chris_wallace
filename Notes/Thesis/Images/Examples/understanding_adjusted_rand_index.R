#!/usr/bin/env Rscript

set.seed(1)
K <- 3
n <- 1e3
x <- c(rep(1, 3874), sample(1:K, size = n, replace = T))
y <- c(rep(1, 3874), sample(1:K, size = n, replace = T))



# comb_data <- data.frame(X = x, Y = y)
# 
# cont_table <- matrix(0, nrow = K, ncol = K)
# for(i in 1:K){
#   for(j in 1:K){
#     cont_table[i, j] <- sum(x == i & y == j) / n
#   }
# }
# 
# cont_table
(cont_table_2 <- table(x, y) / n) 

mcclust::arandi(x , y)
fossil::rand.index(x, y)
fossil:::adj.rand.index(x, y)
