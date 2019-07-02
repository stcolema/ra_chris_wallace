
# Source the function to find gene sets from the KEGG data
source("./Analysis/read_in_kegg_sets.R")

# Read in the universe of genes for this project
probe_key <- fread("./Analysis/probe_key_unique_names.csv")
universe <- probe_key$Gene
my_universe <- probe_key$Unique_gene_name
reduced_universe <- unique(universe)

# The KEGG data
kegg_data <- "Data/kegg_msigdb.txt"

# The number of genes desired in the final gene sample
n_desired <- 625

# The pathways to include
pathway_names <- c(
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_INOSITOL_PHOSPHATE_METABOLISM",
  "KEGG_INFLAMMATORY_BOWEL_DISEASE"
)

# Set a seed for reproducibility
set.seed(1)

# This pathway is not present in the KEGG download so we had to download it 
# separately and parse it with python
missing_path <- pathway_names[3]

# Find the gene sets and reduce to those deisred
gene_sets <- find_gene_sets(kegg_data, universe) #, unique_universe = my_universe)
rel_gene_sets <- gene_sets[names(gene_sets) %in% pathway_names]

# The IBD set is missing from the above, read-in the file created by calling 
# python3 find_ibd_genes.py ./kegg_ibd_genes.txt ibd_genes.csv
ibd_set <- fread("Data/ibd_genes.csv", header = T) %>% unlist() %>% unname()

# Check same genes as expected
# ibd_set[ibd_set %in% universe] %>% sort() 
# my_universe[universe %in% ibd_set]

# Add to our list of gene sets using the unique labels
# rel_gene_sets[[missing_path]] <- my_universe[universe %in% ibd_set]
rel_gene_sets[[missing_path]] <- universe[universe %in% ibd_set]

# Reduce to a vector
genes_to_inclue <- rel_gene_sets %>% unlist() %>% unname() %>% unique()

# The complement to the above set
# other_possibilies <- my_universe[! my_universe %in% genes_to_inclue]
other_possibilies <- universe[! universe %in% genes_to_inclue]

# The number of genes required to get up to n_desired
n_to_sample <- n_desired - length(genes_to_inclue)

# Sample a bunch of random genes
other_genes <- sample(other_possibilies, n_to_sample)

# The final set of genes used
final_genes_used <- c(
  genes_to_inclue,
  other_genes
)

# Change the format to be more in-keeping with previous subsets
newnet <- data.frame(external_gene_name = final_genes_used)

# Save the data in RData format
save(newnet, file = "./Data/1000_gene_set.RData")

