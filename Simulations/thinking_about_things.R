
library(pheatmap)

my_dir <- "~/Desktop/Sim_realistic_test/"
data_file <- "Sim_realistic_k_7_n_1000_p_500.csv"
cluster_data_file <- "Cluster_IDs.csv"
my_data <- read.csv(paste0(my_dir, data_file), row.names = 1)
cluster_IDs <- read.csv(paste0(my_dir, cluster_data_file), row.names = 1)

n <- nrow(my_data)
p <- ncol(my_data)
features_sampled <- floor(sqrt(p))
n_sampled <- floor(n * (1 - 1/exp(1)))

n_seeds <- 10

set.seed(1)

seeds_used <- sample(1:1e6, size = n_seeds)

data_used <- list()
orderings <- list()

for(i in 1:n_seeds){
  set.seed(seeds_used[i])
  .col_used <- sample(1:p, size = features_sampled)
  orderings[[i]] <- .sample_order <- sample(1:n, size = n)
  data_used[[i]] <- .curr_data <- my_data[.sample_order, .col_used]
}
all(row.names(cluster_IDs) == row.names(data_used[[1]])[order(orderings[[1]])])


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

pheatmap(cor(my_data))

plot(cor(my_data) %>% abs() %>%  rowSums())

plotData(data_used[[1]][order(orderings[[1]]), ], cluster_IDs[,1])
plotData(data_used[[2]][order(orderings[[2]]), ], cluster_IDs[,1])
plotData(data_used[[3]][order(orderings[[3]]), ], cluster_IDs[,1])
plotData(data_used[[2]][order(orderings[[2]]), ], cluster_IDs[,1])

plotData(data_used[[2]][order(orderings[[2]]), ], cluster_IDs[,1])

# The results for the real data are not shocking as they have cleaned the data!
real_data <- read.csv("../Data/BCC_data/MDI_data/MDI_Expression_data.csv", row.names = 1)
pheatmap(real_data)
# 
# 
# plot(cor(real_data) %>% abs() %>%  rowSums())
# 
# pheatmap(cor(t(real_data)) %>% abs() %>% rowSums() %>% t(), cluster_cols = T, cluster_rows = F)
# pheatmap(cor(t(real_data)))
# 
# order <- cor(t(real_data)) %>% abs() %>% rowSums() %>% dist() %>%  hclust()
# 
# cor_data <- cor(t(real_data)) %>% abs() %>%  rowSums()
# plot(sort(cor_data))
# plot(cor(t(real_data))[order$order, order$order] %>% abs() %>%  rowSums())

# pheatmap(cor(t(real_data))[order$order, order$order] %>% abs() %>%  rowSums())
