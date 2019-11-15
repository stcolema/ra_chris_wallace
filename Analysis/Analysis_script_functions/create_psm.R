#!/usr/bin/env Rscript

# Functions to create a Posterior similarity matrix

library(data.table, quietly = T)
library(magrittr, quietly = T)
library(Rcpp, quietly = T)
# Rcpp::sourceCpp("/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/posterior_sim_mat.cpp") # install.packages("Rcpp", dep = T)


# find the similarity between the columns of a data.table
psm_entries <- function(my_dt, n, p){
  # n <- nrow(my_dt)
  # p <- ncol(my_dt)
  
  # Declare output list (this will be a lower triangular matrix)
  out <- list()
  
  # Loop over columns
  for(i in 1:(p -  1)){
    out[[i]] <- compare_col_for_psm(i, my_dt)
  }
  
  # Find the proportion of similarity
  lapply(out, `/`, n)
  
}

# Compare the ith column of a data.table to the [i+1,...,n]th columns
compare_col_for_psm <- function(i, my_dt){
  
  # The column of interest
  x <- my_dt[, ..i] %>% c() %>% unlist()
  
  # The indices of the columns ot drop
  drop <- 1:i
  
  # A fill to ensure the output is the same lenght for every column (possibly
  # unnecessary)
  fill <- rep(0, i)
  
  # The number of entries exactly the same between the ith column and the 
  # [i+1,...,n]th columns
  y <- (x == my_dt[, -..drop]) %>% colSums()
  
  # Fill with 0s
  c(fill, y)
}

# Turn a triangular matrix with no diagonal into a symmetric matrix with 1's 
# along the diagonal
construct_symmetrix <- function(tiangular_mat, n, diag_value = 1){
  tiangular_mat + t(tiangular_mat) + (diag(n) * diag_value)
}

# Bind the list of vectors from compare_col_for_psm into a triangular matrix
# Add a column of 0s on the end
construct_psm <- function(list_sim, n = NULL){

  # Bind the entries of the list into a matrix (each entry forms a column)
  psm <- do.call(cbind, list_sim)
  
  if(is.null(n)) {
    n <- ncol(psm)
  }

  # Add an empty column
  psm <- cbind(psm, rep(0, n) )
  psm
  
}

# Wrapper function to find the PSM for a data.table
make_psm <- function(my_dt) {
  
  n <- nrow(my_dt)
  p <- ncol(my_dt)
  
  list_sim <- psm_entries(my_dt, n, p)
  psm <- construct_psm(list_sim, p)
  construct_symmetrix(psm, p)
}

n <- 4
n_iter <- 10
z1 <- matrix(sample(1:12, size = n*n_iter,  replace = T), nrow = n, ncol = n_iter) %>%
  t()
# 
z2 <- z1 %>%
  as.data.table()
# 
# stm_i <- Sys.time()
# 
psm1 <- make_psm(z2)
# stm_j <- Sys.time()
# 
# psm2 <- similarity_mat(t(z1))
# 
# stm_k <- Sys.time()
# 
# if((psm1 - psm2)%>% sum() > 1e-5){
#   print("Error")
# }
# 
# print("Time for data.table psm")
# print(stm_j - stm_i)
# 
# print("Time for C++ psm")
# print(stm_k - stm_j)
# 
# 
# (psm1 %>% as.matrix()) - psm2
#  - psm2
# matrix(psm1) - psm2
# str(psm1)
# str(psm2)
# 
# psm1 %>% unlist() %>% c() %>% matrix(ncol = n) %>% `-`(psm2)
# 
# psm1[1:10, 1:10] - psm2[1:10, 1:10]
# 
# x1 <- matrix(1:4, nrow = 2)
# x2 <- matrix(c(1, 3, 2, -4), nrow = 2)
# x2 - x1
# x2 %>% 
#   `-`(x1)
