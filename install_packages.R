#!/usr/bin/env Rscript

# Script to install the required packages if not already installed

required_packages <- c(
  "tidyverse",
  "mcclust",
  "mclust",
  "pheatmap",
  "Rcpp",
  "RcppArmadillo",
  "lpSolve",
  "data.table",
  "optparse",
  "magrittr",
  "RColorBrewer",
  "png",
  "gridExtra"
)

for (package in required_packages) {
  if (!(package %in% rownames(installed.packages()))) {
    install.packages(package, dep = T)
  }
}
