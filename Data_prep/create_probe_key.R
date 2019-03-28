#!/usr/bin/env Rscript
# File to create a .csv file for the info on the probes present

library(data.table)
# setwd("~/Desktop/MDI/Data")

# Website for probe information
data_url <- paste0(
  "http://139.165.108.18/srv/genmol/permanent/",
  "1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/",
  "CEDAR_GE/Probes_good_reanno_31137_TSS.txt"
)

# Read in probe data from the website
probe_key <- read.csv(data_url, header = T, sep = " ")

# Read in the data of the probes present in our data
probes_present <- fread("Analysis/probe_IDs_present.csv", header = T) %>%
  c() %>%
  unlist() %>%
  unname()

# Take the subset of the probe_key info that is relevant to our data
rel_probes <- probe_key[probe_key$ProbeID %in% probes_present, ] %>%
  .[match(probes_present, .$ProbeID), ]

# Check that the relevant probe information df is in the correct order
if (any(rel_probes$ProbeID != probes_present)) {
  cat("Issue - check order is correct.\n")
  stop("Probe IDs not the same in both sets. Please inspect.")
}

# Visually inspect the heads of both
head(rel_probes)
head(probes_present)

# All is good so let's write to file
fwrite(rel_probes, file = "probe_key.csv")
