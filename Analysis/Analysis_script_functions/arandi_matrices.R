
library(pheatmap)
library(magrittr)
library(mcclust)
library(tibble)


function_dir <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/Analysis_script_functions/"

function_scripts <- c(
  "plot_rand_index.R", # for unlisted arandi
  "create_psm.R" # for constructing symmetric matrices
)

for (f in paste0(function_dir, function_scripts)) {
  source(f)
}

arandi_for_vectorisation <- function(my_vec, break_col, end_col){
  
  my_vec_ <- unlist(my_vec)
  
  arandi(my_vec_[1:break_col], my_vec_[(break_col+1) : end_col])
  
}

arandi_matrices <- function(cl_alloc_list, n_datasets){
  n_iter <- nrow(cl_alloc_list[[1]])
  
  out_list <- list()
  base_matrix <- matrix(0, nrow = n_datasets, ncol = n_datasets)
  
  for(i in 1:n_iter){
    
    for(j in 1 : (n_datasets - 1)) {
      
      for(k in (j + 1) : n_datasets){
        
        cl_j <- cl_alloc_list[[j]][i,] %>% unlist()
        cl_k <- cl_alloc_list[[k]][i,] %>% unlist()
        
        base_matrix[j, k] <- arandi(cl_j, cl_k)
        
      }
    }
    out_list[[i]] <- construct_symmetrix(base_matrix, n_datasets, diag_value = 0)
  }
  out_list
}

average_matrix <- function(matrix_list){
  matrix_list %>% 
    simplify2array() %>% 
    apply(1:2, mean)
}

# === From when we played around ===============================================


if(F){

find_adj_rand_ind <- unlist_arandi

# find_adj_rand_ind  <- function(v1, v2){
#   v1_ <- unlist(v1)
#   v2_ <- unlist(v2)
#   
#   arandi(v1, v2)
# }

x <- readRDS("/home/MINTS/sdc56/Desktop/MDI_small_geneset_outputs/Small_std_rnd_long_18_943700/compare_tibble.rds")
x <- readRDS("/home/MINTS/sdc56/Desktop/MDI_small_geneset_outputs/Small_std_0_1000/compare_tibble.rds")

cl_alloc <- x$mdi_allocation

all_datasets <- unique(x$dataset)

n_datasets <- length(unique(x$dataset))
arandi_list <- list()

t_1 <- Sys.time()

for(i in 1: (n_datasets - 1)){

  arandi_list[[i]] <- list()

  d1 <- cl_alloc[[i]]

  break_col <- ncol(d1)

  for(j in (i + 1) : n_datasets){
    d2 <- cl_alloc[[j]]

    end_col <- ncol(d2) + break_col
    d_lie <- cbind(d1, d2)

    arandi_list[[i]][[j]] <- apply(d_lie, 1, arandi_for_vectorisation, break_col, end_col)
    #
    #
    # a <- find_adj_rand_ind(d1[i,], d2[i,])
    #
  }

}
# str(arandi_list)
# 
# x18 <- do.call(rbind, arandi_list)
# 
# x18[1 ,4][[1]][1]
# x21 <- as.array(x18)
# 
# x21[, ]
# 
# 
# 
# toy_example <- list(
#   list(NULL, c(3, 2, 3), c(-1, 0, 11)),
#   list(NULL, c(11, -2, 3), c(1, 0, 11)),
#   list(NULL, NULL, c(1, 4, 9))
# )
# 
# toy_example %>% str()
# 
# new_toy <- do.call(rbind, toy_example)
# 
# new_toy2 <- new_toy
# 
# new_toy2[] <- lapply(new_toy, )
# 
# lapply(arandi_list[[1]], `[`, 1)
# purrr::map(arandi_list[[1]], 1)
# 
# x17 <- sapply(arandi_list, map, 1) %>% rbind()
# 
# lmap(arandi_list, map, 1)
# 
# unlist(x17)

# t_2 <- Sys.time()
# 
# # n_iter <- nrow(x$mdi_allocation[[1]])
# t_3 <- Sys.time()
# 
# ari_list_2 <- arandi_matrices(cl_alloc, n_datasets)
# # base_matrix <- matrix(0, nrow = n_datasets, ncol = n_datasets)
# 
# # ari_list_2[[1]]
# 
# t_4 <- Sys.time()
# 
# print(t_2 - t_1)
# print(t_4 - t_3)
# 
# # str(ari_list_2)
# # 
# # ari_list_2[[1]] %>% pheatmap()
# # 
# # simplify2array(ari_list_2)
# 
# x34 <- apply(simplify2array(ari_list_2), 1:2, mean) %>% set_colnames(all_datasets) %>% set_rownames(all_datasets)
# x45 <- (Reduce("+", ari_list_2) / length(ari_list_2)) %>% set_colnames(all_datasets)%>% set_rownames(all_datasets)
# 
# diag(x34) <- 0.3
# 
# pheatmap(x34, cluster_rows = F, cluster_cols = F)
# pheatmap(x45)
# 
# str(arandi_list)
# 
# arandi_list[[1]][[2]][1]
# arandi_list[[2]][[3]][1]
# 
}