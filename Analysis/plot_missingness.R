

library(ggplot2)
library(dplyr)
library(magrittr)
library(pheatmap)
library(data.table)
library(tidyr)
library(reshape2)

probes_present_file <- "/home/MINTS/sdc56/Desktop/Full_data_29_05_19/Meta_data/probes_present_per_dataset.csv"
probe_key_file <- "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/probe_key_unique_names.csv"


# All possible datasets (and the names of the probes present columns)
all_datasets <- c(
  "CD14",
  "CD15",
  "CD19",
  "CD4",
  "CD8",
  "IL",
  "PLA",
  "RE",
  "TR"
)

# Read in the file relating the probe IDs to the related gene
probe_key <- fread(probe_key_file, header = T)

probes_present_dt <- fread(probes_present_file, header = T) %>% 
  set_colnames(c("V1", all_datasets))

# Pull out the Gene IDs in the correct order
gene_id <- probe_key$Unique_gene_name
probe_names <- probe_key$ProbeID

probes_present_dt <- probes_present_dt[probes_present_dt$V1 %in% probe_names, ]

# Add the gene names to the probes_present dataframe
probes_present_dt$Gene_names <- gene_id[which(probes_present_dt$V1 %in% probe_names)]

probes_present_dt <- probes_present_dt %>% 
  mutate(All = (CD14 + CD15 + CD19 + CD4 + CD8 + IL + PLA + RE + TR == 9),
         All_not_PLA = (CD14 + CD15 + CD19 + CD4 + CD8 + IL + RE + TR == 8))


gathered_data <- gather(probes_present_dt, Dataset, Present, -V1, -Gene_names)

gathered_data$Dataset <- gathered_data$Dataset %>% as.factor

# gathered_data[gathered_data == T] <- "Present"
# gathered_data[gathered_data == F] <- "Empty"
# 
# gathered_data$value <- gathered_data$value %>% as.factor

summary(gathered_data)

ggplot(gathered_data, aes(Dataset, ..count..)) +
  geom_bar(aes(fill = Present), position = "dodge") +
  labs(title = "Probe presence across datasets",
       y = "Count")
