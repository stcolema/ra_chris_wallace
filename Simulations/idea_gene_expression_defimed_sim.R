

library(MASS)
library(pheatmap)
library(ggplot2)
library(viridis)
library(magrittr)

#' @title Generate dataset
#' @description Generate a dataset based upon the cluster means
#' (assumes each feature is independent)
#' @param cluster_means A k-vector of cluster means defining the k clusters.
#' @param n The number of samples to generate in the entire dataset.
#' @param p The number of columns to generate in the dataset.
#' @param pi A k-vector of the expected proportion of points to be drawn from
#' each distribution.
#' @param row_names The row names of the generated dataset.
#' @param col_names The column names of the generated dataset.
generateDataset <- function(cluster_means, n, p, pi,
                            row_names = paste0("Person_", 1:n),
                            col_names = paste0("Gene_", 1:p)) {
  
  # The number of distirbutions to sample from
  K <- length(cluster_means)
  
  # The membership vector for the n points
  cluster_IDs <- sample(K, n, replace = T, prob = pi)
  
  # The data matrix
  my_data <- matrix(nrow = n, ncol = p)
  
  # Iterate over each of the columns permuting the means associated with each
  # label.
  for (j in 1:p)
  {
    reordered_cluster_means <- sample(cluster_means)
    
    # Draw n points from the K univariate Gaussians defined by the permuted means.
    for (i in 1:n) {
      my_data[i, j] <- rnorm(1, mean = reordered_cluster_means[cluster_IDs[i]])
    }
  }
  
  # Order based upon allocation label
  row_order <- order(cluster_IDs)
  
  # Assign rownames and column names
  rownames(my_data) <- row_names
  colnames(my_data) <- col_names
  
  # Return the data and the allocation labels
  list(
    data = my_data[row_order, ],
    cluster_IDs = cluster_IDs[row_order]
  )
}

# This is an idea for generating data
generateFullDataset <- function(K, n, p, p_noisy = 0, a = 2, b = 2) {
  cluster_means <- (1:K - ceiling(K / 2))
  pi <- rbeta(K, a, b)
  
  
  my_data <- generateDataset(cluster_means, n, p, pi)
  
  my_data
}

# Badly named heatmapping function
plotData <- function(x, cluster_IDs,
                     col_pal = colorRampPalette(c("#146EB4", "white", "#FF9900"))(100),
                     my_breaks = mdiHelpR::defineDataBreaks(x, col_pal, mid_point = 0),
                     main = "gen_dataset",
                     ...) {
  
  anno_row <- data.frame(Cluster = factor(paste("Cluster", cluster_IDs))) %>%
    set_rownames(rownames(x))
  
  K <- length(unique(cluster_IDs))
  
  ann_colours <- list(Cluster = viridis(K) %>%
                        set_names(paste("Cluster", sort(unique(cluster_IDs)))))
  
  pheatmap(x,
           color = col_pal,
           breaks = my_breaks,
           annotation_row = anno_row,
           annotation_colors = ann_colours,
           main = main,
           ...
  )
}

genDataFromTable <- function(my_table){
  data_lst <- list()
  for(i in 1:nrow(my_table)){
    K <- my_table[i,]$K
    
    
    n <- my_table[i,]$n
    p <- my_table[i,]$p
    p_signal <- max(0.2 * p, min(p, 100))
    p_noisy <- p - p_signal
    
    a <- 0.5
    b <- 0.5
    
    if(my_table[i,]$mu == "1:K"){
      cluster_means <- 1:K - ceiling(K / 2)
    } else {
      cluster_means <- rnorm(K)
    }
    
    if(my_table[i,]$pi == "constant"){
      pi <- rep(1 / K, K)
    } else {
      pi <- rbeta(K, a, b)
    }
    
    data_lst[[i]] <-  generateDataset(cluster_means, n, p_signal, pi)
    
    .curr_data <- data_lst[[i]]$data
    
    if(p_noisy > 0){
      min_data <- min(cluster_means)
      max_data <- max(cluster_means)
      
      noisy_data <- data.frame(
        matrix(
          rnorm(n * p_noisy, mean = mean(cluster_means), sd = sd(.curr_data)),
          ncol = p_noisy
        )
      ) %>%
        set_rownames(row.names(.curr_data)) %>%
        set_colnames(paste0("Noise_", 1:p_noisy))
      
      data_lst[[i]]$data <- cbind(.curr_data, noisy_data)
    }
  }
  data_lst
}

# === Simulations ==============================================================

# define our ggplot2 theme of choice
theme_set(theme_bw() +
            theme(strip.background = element_rect(fill = "#21677e")) +
            theme(strip.text = element_text(colour = "white")))

# Set a seed for reproducibility
set.seed(1)

# Our colour palette of choice of heatmaps
col_pal <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(100)

# Breaks for correlation plots
cor_breaks <- mdiHelpR::defineBreaks(col_pal)


n <- 100
n1 <- n * 0.3
n2 <- n * 0.1
n3 <- n * 0.1
n4 <- n * 0.5

p_1 <- 1
p_2 <- 2
p_3 <- 3
p_4 <- 2
p_5 <- 2
p_6 <- 1

p <- c(
  p_1,  p_2,  p_3,  p_4,  p_5, p_6
)


my_data <- runif(n, min = -2, max = 2)


for(i in 1:p_2){
  my_data <- cbind(my_data, c(rnorm(n1, mean = 1), rnorm(n - n1, mean = -1)))
}

for(i in 1:p_3){
  my_data <- cbind(my_data, c(rnorm(n1, mean = -2), rnorm(n2 + n3, mean = 2), rnorm(n4, mean = 0)))
}

for(i in 1:p_4){
  my_data <- cbind(my_data, c(rnorm(n1 + n2, mean = 1), rnorm(n3 + n4, mean = -1)))
}

for(i in 1:p_5){
  my_data <- cbind(my_data, c(rnorm(n1, mean = 0),
                              rnorm(n2, mean = 2),
                              rnorm(n3 + n4, mean = -2)
  )
  )
}

for(i in 1:p_6){
  my_data <- cbind(my_data, c(rnorm(n1, mean = -1),
                              rnorm(n2, mean = 0),
                              rnorm(n3, mean = 2),
                              rnorm(n4, mean = -2)
  )
  )
}

my_data <- my_data %>% 
  set_colnames(paste0("Gene_", 1:sum(p))) %>% 
  set_rownames(paste0("Person_", 1:n)) %>% 
  as.data.frame()

# my_data <- matrix(
#   c(
#     col_1,
#     col_2,
#     col_3,
#     col_4,
#     col_5
#   ), ncol = 5
# ) 

annotation_row <- data.frame(Cluster = as.factor(c(rep(1, n1), rep(2, n2), rep(3, n3), rep(4, n4)))) %>% 
  set_rownames(paste0("Person_", 1:n))

annotation_colours <- list(Cluster = 
                             viridis(4) %>% set_names(1:4))


my_breaks <- mdiHelpR::defineDataBreaks(my_data, col_pal)

pheatmap(my_data, cluster_cols = F,cluster_rows = F,
         annotation_row = annotation_row,
         annotation_colors = annotation_colours,
         color = col_pal,
         breaks = my_breaks)

my_data %>% 
  t() %>% 
  cor() %>% 
  pheatmap(annotation_row = annotation_row,
           annotation_colors = annotation_colours,
           color = col_pal,
           breaks = cor_breaks)
