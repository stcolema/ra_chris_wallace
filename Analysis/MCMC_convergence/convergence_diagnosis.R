#!/usr/bin/env Rscript

# === Libraries ================================================================

# devtools::install_github("sbfnk/fitR")
# devtools::install_github("thomasp85/patchwork")

library(magrittr)
library(rjags)
library(fitR) # some extensions for rjags
library(ggplot2)
library(coda)
library(patchwork)

# === Functions ================================================================

# Given a string of the format "aaaa_111" extract the numbers (in this case "111")
extractNumericalFromString <- function(my_string) {
  gsub("[^0-9.]", "", my_string)
}

# Calculate the expected sample size for a chain for a given test burn-in
findESS <- function(mcmc_out, test_burn_in) {

  # initialise data.frame of ess estimates
  ess_burn_in <- data.frame(t(effectiveSize(mcmc_out)))

  # loop over all test values after 0
  for (burn_in in test_burn_in[-1]) {

    # test burn-in
    test.trace <- burnAndThin(mcmc_out, burn = burn_in)

    # estimate ESS and at to vector of ess estimates
    ess_burn_in <- rbind(ess_burn_in, t(effectiveSize(as.mcmc(test.trace))))
  }

  ess_burn_in$burn_in <- test_burn_in

  ess_burn_in
}

# From Stack Exchange https://stackoverflow.com/questions/22908050/formula-evaluation-with-mutate
# Use dplyr::mutate with a formula encoded as a string (used as we do not want
# to restrict the number of variables present)
stringMutate <- function(df, str_formula) {
  q <- quote(mutate(df, Total = str_formula))
  eval(parse(text = sub("str_formula", str_formula, deparse(q))))
}

# Create parameter names for output from MDI for a given number of datasets
createParameterNames <- function(n_datasets) {
  # The parameters of interest

  # The mass parameters - 1 per dataset
  mass_param_names <- paste0("MassParameter_", 1:n_datasets)

  # The phi parameters are combinatorical
  phi_indices <- combn(1:n_datasets, 2) %>%
    apply(., 2, paste0, collapse = "")

  phi_param_names <- paste0("Phi_", phi_indices)

  all_parameters <- c(mass_param_names, phi_param_names)

  all_parameters
}

# Find the appropriate burn-in for a given chain of MCMC (approximate and
# possible to fail - check plottd output to be sure)
findBurnIn <- function(mcmc_out,
                       n_datasets,
                       max_burn_in = ifelse(is.data.frame(mcmc_out) | is.mcmc(mcmc_out), nrow(mcmc_out), nrow(mcmc_out[[1]])) / 2,
                       step_size = round(max_burn_in / 50)) {


  # What parameters are present
  all_parameters <- createParameterNames(n_datasets)

  # The function describing the sum of all parameters ESS
  param_func <- paste(all_parameters, collapse = "+")


  # Vector of values to test effective sample size (ESS) at
  test_burn_in <- seq(0, max_burn_in, step_size) # test values

  # If not already a list, change mcmc_out to a list
  if (!class(mcmc_out) %in% c("mcmc.list", "list")) {
    mcmc_out <- list("chain_1" = mcmc_out)
  }

  # If there are no names, set these to chain_{1..n_chains}
  if (is.null(names(mcmc_out))) {
    names(mcmc_out) <- paste0("chain_", seq_along(mcmc_out))
  }

  # Create a data.frame of effective sample size for each parameter at each
  # possible burn in
  df_ess.burn.in <- plyr::ldply(mcmc_out, .fun = findESS, test_burn_in)

  # Create a variable comparing the ESS for each variable
  est_burn_in_df <- df_ess.burn.in %>%
    stringMutate(param_func)

  # Choose the burn in for which the parameters have the largest combined ESS
  est_burn_in <- est_burn_in_df$burn_in[est_burn_in_df$Total == max(est_burn_in_df$Total)]

  est_burn_in
}

# Save ESS plots
saveESSPlots <- function(my_data, n_seeds, save_dirs, plot_titles,
                         gen_ess_title = "Burn in diagnosis plot",
                         plot_type = ".png") {
  ess_file_names <- paste0(save_dirs, "/burn_in_plot", plot_type)
  ess_titles <- paste0(plot_titles, ": ", gen_ess_title)

  for (i in 1:n_seeds) {

    # Plot title and filename
    curr_file_name <- ess_file_names[i]
    curr_title <- ess_titles[i]

    # Create plot
    plotESSBurn(my_data[[i]]) +
      labs(
        title = curr_title,
        x = "Tested burn in",
        y = "Effective sample size"
      )

    # Save
    ggsave(curr_file_name)
  }
}

# Save autocorrelation plots using ggplot2 style
saveAutoCorrelationPlots <- function(mcmc_data, n_seeds, save_dirs, plot_titles,
                                     lag.max = NULL,
                                     gen_auto_corr_title = "Autocorrelation",
                                     plot_type = ".png") {

  # The filenames and plot titles for the plots
  auto_corr_file_names <- paste0(save_dirs, "/autocorrelation_plot", plot_type)
  auto_corr_titles <- paste0(plot_titles, ": ", gen_auto_corr_title)

  for (ii in 1:n_seeds) {
    x <- mcmc.list(as.mcmc(mcmc_data[[ii]]))

    # Plot title and filename
    curr_file_name <- auto_corr_file_names[ii]
    curr_title <- auto_corr_titles[ii]


    # Consider the chain as a time series and calculate the auto-correlation
    x_auto_cor <- x[[1]] %>%
      as.ts() %>%
      acf(lag.max = lag.max, plot = FALSE)

    # Create a data.frame to hold the autocorrelation info
    auto_cor_df <- data.frame(Lag = x_auto_cor$lag[, 1, 1])

    # Add the autocorrelation for each variable to said data.frame
    for (jj in 1:nvar(x)) {
      auto_cor_df <- cbind(auto_cor_df, x_auto_cor$acf[, jj, jj])
    }

    colnames(auto_cor_df) <- c("Lag", varnames(x))

    # The variables to be considered in going from wide to long
    vars_to_gather <- varnames(x)

    # The long form of the original dataset (to better fit ggplot2)
    auto_cor_df_long <- tidyr::gather_(auto_cor_df,
      "Parameter",
      "Autocorrelation",
      vars_to_gather,
      factor_key = TRUE
    )

    # Plots!
    p_lst[[ii]] <- ggplot(
      auto_cor_df_long,
      aes(x = Lag, y = Autocorrelation, colour = Parameter)
    ) +
      geom_line() +
      labs(
        title = paste0("Seed ", ii, ": ", gen_auto_corr_title)
      )

    ggsave(curr_file_name)
  }
  # Return a list of plots
  p_lst
}

# === Setup ====================================================================

# Set ggplot2 theme
theme_set(theme_bw())

# Example directory for output data
gen_mdi_dir <- "~/Documents/PhD/Year_1/Consensus_clustering/Yeast/Output_data/MDI/"

# List the sub directories
sub_dirs <- list.dirs(
  path = gen_mdi_dir,
  full.names = F,
  recursive = FALSE
)

# The generic file name before the seed and extension
gen_file_name <- "out_"

# The file output type (this should never change)
file_ext <- ".csv"

# One of .png or .pdf
plot_type <- ".png"

# Instructions regarding plots to construct
plot_burn_in <- TRUE
plot_auto_corr <- TRUE
plot_geweke <- TRUE
plot_gelman <- TRUE

# The file names will be in order of strings - 10 comes before 2. Remedy this.
numerical_order <- sub_dirs %>%
  str_extract_numerical() %>%
  as.numeric() %>%
  order()

# Rearrange the sub directories to be intuitive
sub_dirs <- sub_dirs[numerical_order]
full_sub_dirs <- paste0(gen_mdi_dir, sub_dirs)

# This is the sub directory we will save our diagnostic plots to
new_sub_dirs <- "Convergence_diagnostics"
new_full_sub_dirs <- paste(full_sub_dirs, new_sub_dirs, sep = "/")

# Create these directories
x <- new_full_sub_dirs %>%
  sapply(., dir.create, showWarnings = FALSE)

# The files to be read in
file_names <- paste0(gen_file_name, sub_dirs, file_ext)

# The filename combined iwth the directory they reside in
full_files_names <- paste0(gen_mdi_dir, sub_dirs, "/", file_names)

# The sub directory names converted to sentence case for plot title
pretty_sub_dir_names <- sub_dirs %>%
  stringr::str_to_sentence(locale = "en") %>%
  stringr::str_replace("_", " ")

# The number of chains being compared
n_seeds <- length(sub_dirs)

# the number of datasets used by MDI
n_datasets <- 3

# the number of continuous variables outputted
n_cols_cont <- n_datasets + choose(n_datasets, 2)

# === Read in data =============================================================

# Load the relevant data
my_data <- vector("list", length = n_seeds)

for (i in 1:n_seeds) {

  # The full file name to load
  f <- full_files_names[i]

  # Read in the data for the continuous variables and convert to an MCMC object
  my_data[[i]] <- read.csv(f)[, 1:n_cols_cont] %>%
    mcmc()
}

# === Plotting =================================================================

# If instructed to plot burn in diagnosis
if (plot_burn_in) {
  saveESSPlots(my_data, n_seeds, new_full_sub_dirs, pretty_sub_dir_names,
    gen_ess_title = "Burn in diagnosis plot",
    plot_type = ".png"
  )
}


# Find the common burn in across chains
common_burn_in <- sapply(my_data, findBurnIn, 3) %>%
  max()

# Apply the burn in to the data
new_data <- lapply(my_data, burnAndThin, burn = common_burn_in)

# Plot autocorrelation for each seed for each parameter
if (plot_auto_corr) {
  p_lst <- saveAutoCorrelationPlots(new_data,
    n_seeds,
    new_full_sub_dirs,
    pretty_sub_dir_names,
    lag.max = NULL,
    gen_auto_corr_title = "Autocorrelation",
    plot_type = plot_type
  )
}

p_lst[[1]] / p_lst[[10]] + plot_layout(guides = "collect")

gelman.plot(new_data)
