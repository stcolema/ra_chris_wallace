

library(data.table)
library(magrittr)

cont_to_cat <- function(x, bounds){
  
  # Assign the value to a given category depending on the percentile
  if(x < bounds[1]){
    x <- "Zero"
    return(x)
  }
  
  if(x < bounds[2]){
    x <- "Small"
    return(x)
  }
  if(x < bounds[3]){
    x <- "Medium"
    return(x)
  }
  x <- "Large"
  x
}

# File path to data files
file_path <- "/home/MINTS/sdc56/Desktop/Matlab_input_small_names"

# The file names
file_names <- list.files(path = file_path, full.names = T, include.dirs = F) %>%
  grep("csv", ., value = TRUE)

# Stripped of file path and extension
cleaned_file_names <- file_names %>% 
  tools::file_path_sans_ext() %>% 
  basename()

n_datasets <- length(cleaned_file_names)

# The output list
cat_data <- list()

# The names of each dataset
dataset_names <- cleaned_file_names %>%
  gsub("([^\\_]+)\\_.*", "\\1", .)

generic_output_name <- paste0(file_path, "/categorical_data_")

# Loop over the datasets saving a new version of each one with categorical variables
for(j in 1:n_datasets){
  
  curr_dataset <- fread(file_names[[j]])
  
  na_curr_dataset <- curr_dataset
  
  na_curr_dataset[na_curr_dataset == 0] <- NA
  
  bounds <- apply(na_curr_dataset[,-1], 2, quantile, c(0, 1/3, 2/3, 3/3), na.rm = T)
  
  for(i in 2:ncol(curr_dataset)){
    curr_dataset[,i] <- apply(curr_dataset[,..i], 1, cont_to_cat, bounds[,(i-1)]) %>% 
      as.factor()
  }
  
  curr_name <- dataset_names[[j]]
  
  curr_output_name <- paste0(generic_output_name, curr_name, ".csv")
  
  fwrite(curr_dataset, curr_output_name)
  
  cat_data[[curr_name]] <- curr_dataset
  
}
