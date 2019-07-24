

library(ggplot2)
library(dplyr)
library(magrittr)
library(pheatmap)
library(data.table)
library(tidyr)
library(reshape2)

add_level <- function(x, new_level){
  if(is.factor(x)) return(factor(x, levels=c(levels(x), new_level)))
  return(x)
}

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

# probes_present_dt <- probes_present_dt %>% 
#   mutate(All = (CD14 + CD15 + CD19 + CD4 + CD8 + IL + PLA + RE + TR == 9),
#          `All (excl. PLA)` = (CD14 + CD15 + CD19 + CD4 + CD8 + IL + RE + TR == 8))


gathered_data <- gather(probes_present_dt, Dataset, Present, -V1, -Gene_names)

gathered_data$Dataset <- gathered_data$Dataset %>% as.factor

# gathered_data[gathered_data == T] <- "Present"
# gathered_data[gathered_data == F] <- "Empty"
# 
# gathered_data$value <- gathered_data$value %>% as.factor

summary(gathered_data)
# gathered_data$Dataset <- add_level(gathered_data$Dataset, "All datasets excl. PLA")
# gathered_data$Dataset[is.na(gathered_data$Dataset)] <- "All datasets excl. PLA"


new_order <- levels(gathered_data$Dataset)[c(1,3:6, 8:9, 2, 7)] #:8, 10:11, 9, 1, 2)]

gathered_data <- within(gathered_data,
                   Dataset <- factor(Dataset,
                                      levels=new_order)
)

sum(gathered_data$Present[gathered_data$Dataset == "PLA"])

n_pres <- c()
for(d in new_order[1:9]){
  n_pres <- c(n_pres, sum(gathered_data$Present[gathered_data$Dataset == d]))
}

mean(n_pres)
mean(n_pres[-9])

ggplot(gathered_data, aes(Dataset, ..count.., fill = Present)) +
  geom_bar(, position = "dodge") +
  ylim(0, 20000) +
  geom_text(stat='count', aes(label=..count..), position = position_dodge(width=1), hjust = 0, vjust=-1, angle = 60) +
  labs(title = "Probe presence across datasets",
       y = "Count") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust=1)) 

ggsave("./Notes/Thesis/Images/probe_presence_across_datasets_no_all.png")
