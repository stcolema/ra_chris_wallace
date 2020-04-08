

library(MASS)
library(pheatmap)
library(ggplot2)
library(viridis)
library(magrittr)
library(mdiHelpR)
library(ggfortify)

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


gaussian_means <- matrix(
  c(-3, -3, -3, 3, 3, -3, 3, 3, 0, 0),
  nrow = 2
)

k <- ncol(gaussian_means)
n <- 500
p <- nrow(gaussian_means)

Sigma <- diag(1, nrow = p)

for(i in 1:k){
  
  cluster_data <- mvrnorm(n / k, gaussian_means[, i], Sigma)
  
  if(i == 1){
    my_data <- cluster_data
  } else {
    my_data <- rbind(my_data, cluster_data)
  }
  
}

my_data <- data.frame(
  Gene_1 = my_data[, 1],
  Gene_2 = my_data[, 2]
) %>% set_rownames(paste0("Person_", 1:n))

plot_data <- my_data
plot_data$Cluster <- as.factor(rep(1:k, each = n / k))

ggplot(plot_data, aes(x = Gene_1, y = Gene_2, colour = Cluster)) +
  geom_point() +
  labs(
    title = "Simple mixture of Gaussians",
    x = "Gene 1",
    y = "Gene 2"
  )  +
  scale_colour_viridis_d()


plotData(my_data, rep(1:k, each = n / k),
         main = "Expression data across two genes")


row_order <- order(my_data[, 1])
plotData(my_data[row_order, ], rep(1:k, each = n / k)[row_order],
         main = "Expression data across two genes",
         cluster_rows = F)


# Table of scenarios to generate
my_table <- data.frame(
  n = rep(c(1e2, 1e3), each = 8),
  p = rep(c(30, 500), each = 2, 4),
  K = c(rep(c(3, 5), 4), rep(c(5, 7), 4)),
  pi = rep(c("constant", "varying"), each = 4, 2),
  mu = rep(c("1:K", "rnorm(K)"), 2, each = 2),
  sigma = "I"
)

my_table

# Generate the scenarios captured in our table
data_lst <- genDataFromTable(my_table)

pc_1 <- prcomp(data_lst[[1]]$data)
autoplot(pc_1, data = data_lst[[1]]$data) +
  geom_point(aes(colour = as.factor(data_lst[[1]]$cluster_IDs))) +
  labs(title = "PCA of generated data",
       subtitle = "Coloured by cluster IDs",
       colour = "Cluster") + 
  scale_color_viridis_d()

# Look at these!
plotData(data_lst[[1]]$data, data_lst[[1]]$cluster_IDs, cluster_rows = T)
plotData(data_lst[[2]]$data, data_lst[[2]]$cluster_IDs, cluster_rows = T)
plotData(data_lst[[3]]$data, data_lst[[3]]$cluster_IDs, cluster_rows = T)
plotData(data_lst[[4]]$data, data_lst[[4]]$cluster_IDs, cluster_rows = T)
plotData(data_lst[[5]]$data, data_lst[[5]]$cluster_IDs, cluster_rows = T)
plotData(data_lst[[6]]$data, data_lst[[6]]$cluster_IDs, cluster_rows = T)
plotData(data_lst[[7]]$data, data_lst[[7]]$cluster_IDs, cluster_rows = T)
plotData(data_lst[[8]]$data, data_lst[[8]]$cluster_IDs, cluster_rows = T)
plotData(data_lst[[9]]$data, data_lst[[9]]$cluster_IDs, cluster_rows = T)
plotData(data_lst[[10]]$data, data_lst[[10]]$cluster_IDs, cluster_rows = T)
plotData(data_lst[[11]]$data, data_lst[[11]]$cluster_IDs, cluster_rows = T, cluster_cols = F)
plotData(data_lst[[12]]$data[, 1:100], data_lst[[12]]$cluster_IDs, cluster_rows = T, cluster_cols = T)
plotData(data_lst[[14]]$data, data_lst[[14]]$cluster_IDs, cluster_rows = F, cluster_cols = T)
plotData(data_lst[[15]]$data[,1:100], data_lst[[15]]$cluster_IDs, cluster_rows = T, cluster_cols = T)
plotData(data_lst[[16]]$data[,1:100], data_lst[[16]]$cluster_IDs, cluster_rows = T, cluster_cols = T)

plotData(data_lst[[3]]$data[, 1:100], data_lst[[3]]$cluster_IDs,
         cluster_rows = T)

plotData(data_lst[[4]]$data[, 1:100], data_lst[[4]]$cluster_IDs,
         cluster_rows = T)

# plotData(data_lst[[7]]$data[, 1:100], data_lst[[7]]$cluster_IDs, cluster_rows = T)
plotData(data_lst[[8]]$data[, 1:100], data_lst[[8]]$cluster_IDs, 
         cluster_rows = T)

plotData(data_lst[[11]]$data[, 1:100], data_lst[[11]]$cluster_IDs, cluster_rows = T, cluster_cols = T)
plotData(data_lst[[12]]$data[, 1:100], data_lst[[12]]$cluster_IDs, cluster_rows = T, cluster_cols = T)
plotData(data_lst[[15]]$data[,1:100], data_lst[[15]]$cluster_IDs, cluster_rows = T, cluster_cols = T)
plotData(data_lst[[16]]$data[,1:100], data_lst[[16]]$cluster_IDs, cluster_rows = T, cluster_cols = T)

# data_lst[[16]]$cluster_IDs %>%
#   as.data.frame() %>% 
#   set_rownames(row.names(data_lst[[16]]$data)) %>% 
#   set_colnames("Cluster") %>% 
#   write.csv("~/Desktop/Sim_realistic_test/Cluster_IDs.csv")

data_lst[[11]]$data[, 1:100] %>% 
  t() %>% 
  cor() %>% 
  pheatmap(color = col_pal, breaks = cor_breaks)

data_lst[[12]]$data[, 1:100] %>% 
  t() %>% 
  cor() %>% 
  pheatmap(color = col_pal, breaks = cor_breaks)

data_lst[[15]]$data[,1:100] %>% 
  t() %>% 
  cor() %>% 
  pheatmap(color = col_pal, breaks = cor_breaks)

data_lst[[15]]$data %>% 
  t() %>% 
  cor() %>% 
  pheatmap(color = col_pal, breaks = cor_breaks)


data_lst[[16]]$data[,1:100] %>% 
  t() %>% 
  cor() %>% 
  pheatmap(color = col_pal, breaks = cor_breaks)

data_lst[[16]]$data %>% 
  t() %>% 
  cor() %>% 
  pheatmap(color = col_pal, breaks = cor_breaks)


n <- 100
n1 <- n * 0.3
n2 <- n * 0.1
n3 <- n * 0.1
n4 <- n * 0.5

p_1 <- 1
p_2 <- 2
p_3 <- 2
p_4 <- 2
p_5 <- 2

p <- c(
  p_1,  p_2,  p_3,  p_4,  p_5
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
                              rnorm(n2, mean = rnorm(1, mean = 2)),
                              rnorm(n3 + n4, mean = rnorm(1, mean = -2))
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

pheatmap(my_data, cluster_cols = F,cluster_rows = F,
         annotation_row = annotation_row,
         annotation_colors = annotation_colours,
         color = col_pal)

my_breaks <- mdiHelpR::defineBreaks(col_pal)

my_data %>% 
  t() %>% 
  cor() %>% 
  pheatmap(annotation_row = annotation_row,
           annotation_colors = annotation_colours,
           color = col_pal,
           breaks = my_breaks)
