

library(data.table)

yeast_1 <- fread("/home/MINTS/sdc56/Desktop/Yeast_data/MDItestdata1.csv", header = T)

yeast_genes <- yeast_1[,1]

fwrite(yeast_genes, "/home/MINTS/sdc56/Desktop/Yeast_data/Yeast_genes.csv")

genes_present <- data.table(yeast_genes, D1 = T, D2 = T, D3 = T, D4 = T, D5 = T, D6 = T)
genes_present$Gene_names

yeast_probe_key <- data.table(ProbeID = unlist(yeast_genes), Unique_gene_name = unlist(yeast_genes))

fwrite(yeast_probe_key, "/home/MINTS/sdc56/Desktop/Yeast_data/Yeast_probe_key.csv")

fwrite(genes_present, "/home/MINTS/sdc56/Desktop/Yeast_data/Yeast_genes_present_per_dataset.csv")
