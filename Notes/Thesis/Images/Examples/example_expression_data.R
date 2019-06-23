
library(pheatmap)
library(magrittr)

save_dir <- "/Users/steph/Desktop/Bioinformatics/GEN80436 - Thesis/ra_chris_wallace/Notes/Thesis/Images/Examples/"
dir.create(save_dir, showWarnings = F)

x <- data.frame(`Person 1` = c(5.1, 5.1, 1.4, 1.4, 1.4),
           `Person 2` = c(5.2, 4.9, 1.5, 1.2, 1.5),
           `Person 3` = c(4.9, 5.2, 1.2, 1.5, 1.4),
           `Person 4` = c(5.0, 5.4, 1.3, 1.7, 1.5)
           ) %>% 
  set_rownames(paste0("Gene ", LETTERS[1:5]))

file_1 <- "example_expression_data.png"
file_2 <- "example_standardised_expression_data.png"

pheatmap(x, filename = paste0(save_dir, file_1))

y <- x %>% t() %>% scale() %>% t()

y %>% round(2)

pheatmap(y, filename = paste0(save_dir, file_2))
