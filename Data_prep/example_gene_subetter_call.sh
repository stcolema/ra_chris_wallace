#!/usr/bin/env bash

Rscript Data_prep/gene_subsetter.R -f ../MDI/Data/Fill_NAs_Min/CD4.csv -g small -p Analysis/probe_key.csv --gene_data stephen-genesets.RData -t TRUE

