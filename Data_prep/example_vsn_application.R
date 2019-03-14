#!/usr/bin/env Rscript

# For the pipe
library(magrittr)

# For fread and fwrite
library(data.table)

# For variance stabilisation functions
library(vsn)


vsn_data_table <- function(dt, exponent = 2, id_col = 1){
  
  # Add row names from the id column
  row.names(dt) <- dt[, id_col]
  
  # Remove the id column, remove any empty observations
  # Convert to matrix, return to base values (as log2)
  # Apply vsn
  dt %>%
    .[, -id_col] %>%             # Remove the id column
    .[rowSums(.) != 0,] %>%      # Remove any empty observations
    as.matrix %>%                # Convert to matrix (required by vsnMatrix)
    exponent^. %>%               # Exponentiate to move from log scale
    vsnMatrix                    # Apply variance stabilization
  
}

dt <- read.table("/home/MINTS/sdc56/Desktop/MDI/Data/Fill_NAs_Min/CD4.csv", 
                 header = T, 
                 sep = ",")

# Add row names from the id column
# row.names(dt) <- dt[, id_col]
# 
# row.names(dt) <- dt$V1

my_vsn_data <- vsn_data_table(dt)

meanSdPlot(my_vsn_data, ranks = T)
meanSdPlot(my_vsn_data, ranks = F)

vsn_data <- dt %>% 
  .[,-1] %>% 
  .[rowSums(.) != 0,] %>% 
  as.matrix %>%
  2^. %>%
  vsnMatrix

meanSdPlot(vsn_data, ranks = T)
meanSdPlot(vsn_data, ranks = F)

x2 <- data.table(exprs(my_vsn_data))
