
# To remove the V1 columm from gene subset data
library(magrittr)
library(data.table)

home_dir <- "~/Desktop/ra_chris_wallace/Data/"
dirs_to_read <- c("Big_gene_set", "Med_gene_set", "Small_gene_set") %>% 
  paste0(home_dir, .)

for(curr_dir in dirs_to_read){
  files_present <- list.files(path = curr_dir) %>% 
    grep("vsn_*", ., value = T) %>% 
    paste(curr_dir, ., sep = "/")
  
  for(f in files_present){
    curr_dt <- fread(f, header = T) %>% 
      .[, -1]
    
    fwrite(curr_dt, file = f)
    
  }
}
