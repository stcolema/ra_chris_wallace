
library(ggplot2)
library(pheatmap)
x1 <- read.csv("../Data/MDI_test_data/GaussianTestData1.csv", row.names = 1)
pheatmap(x1)

x2 <- read.csv("../Data/MDI_test_data/MDItestdata1.csv", row.names = 1)
pheatmap(x2)

summary(x2)

pheatmap(scale(x2))
summary(scale(x2))

for(i in 1:6){
  x <- read.csv(paste0("../Data/MDI_test_data/MDItestdata", i, ".csv"), row.names = 1)
  y <- scale(x)
  write.csv(y, file = paste0("../Data/MDI_test_data/ScaledMDItestdata", i, ".csv"))
}


