
set.seed(1)
K <- 3
n <- 1e6
x <- sample(1:K, size = n, replace = T)
y <- sample(1:K, size = n, replace = T)

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
