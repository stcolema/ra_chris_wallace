#!/usr/bin/env bash

# Example call to analysis_script.R
Rscript ra_chris_wallace/Analysis/analysis_script.R -d ./MDI_small_geneset_outputs/matlab_output_no_vsn_no_norm/ --n_iter 445 --datasets "CD14.csv CD19.csv CD4.csv CD8.csv IL.csv RE.csv TR.csv" -b 100 --time TRUE --plot_trees FALSE --plot_rand_index FALSE --plot_heatmap_clusterings FALSE --plot_phi_series FALSE --plot_phi_densities FALSE --plot_type .pdf
