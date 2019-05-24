

# Create an additional column in the probe_key containing unique gene identifiers

library(magrittr)
library(data.table)

# For creating unique gene names
unique_names_from_recurring <- function(name_vec) {
  
  # Number of entries
  n <- length(name_vec)
  
  # Duplicate the name vec - this will hold the transformed versions
  duplicate_vec <- name_vec
  
  for (i in 1:n) {
    n_occurences <- sum(name_vec == name_vec[i])
    
    if (n_occurences > 1) {
      duplicate_vec[i] <- paste0(
        duplicate_vec[i],
        ".",
        sum(name_vec[1:i] == name_vec[i])
      )
    }
  }
  duplicate_vec
}

probe_key_file <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/probe_key.csv"

probe_key <- data.table::fread(probe_key_file)

unique_gene_names <- probe_key$Gene %>% 
  unique_names_from_recurring()

probe_key$Unique_gene_name <- unique_gene_names

fwrite(probe_key, "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/probe_key_unique_names.csv")

