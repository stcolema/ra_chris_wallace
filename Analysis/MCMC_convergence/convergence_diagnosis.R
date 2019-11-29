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
library(dplyr)
library(plyr)
library(tidyr)
library(stringr)

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
  q <- quote(dplyr::mutate(df, Total = str_formula))
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
saveESSPlots <- function(my_data, save_dirs, plot_titles,
                         n_seeds = length(my_data),
                         gen_ess_title = "Burn in diagnosis plot",
                         plot_type = ".png",
                         save_plots = FALSE) {
  ess_file_names <- paste0(save_dirs, "/burn_in_plot", plot_type)
  ess_titles <- paste0(plot_titles, ": ", gen_ess_title)

  p_lst <- vector("list", n_seeds)

  for (i in 1:n_seeds) {

    # Plot title and filename
    curr_file_name <- ess_file_names[i]
    curr_title <- ess_titles[i]

    # Create plot
    p_lst[[i]] <- plotESSBurn(my_data[[i]]) +
      labs(
        title = curr_title,
        x = "Tested burn in",
        y = "Effective sample size"
      )

    # Save
    if (save_plots) {
      ggsave(curr_file_name)
    }
  }
  p_lst
}

autoCorrelationPlot <- function(x, plot_title, lag.max = NULL, facet = FALSE) {

  # The preferred object type for interacting with coda functions
  x <- mcmc.list(x)

  # Consider the chain as a time series and calculate the auto-correlation
  x_auto_cor <- x[[1]] %>%
    as.ts() %>%
    acf(lag.max = lag.max, plot = FALSE)

  # Create a data.frame to hold the autocorrelation info
  auto_cor_df <- data.frame(Lag = x_auto_cor$lag[, 1, 1])

  # Add the autocorrelation for each variable to said data.frame
  for (i in 1:nvar(x)) {
    auto_cor_df <- cbind(auto_cor_df, x_auto_cor$acf[, i, i])
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

  # Create a line plot of the autocorrelation in each parameter
  # Either split parameter's by colour or facet wrap
  if (facet) {
    # Plots!
    p <- ggplot(
      auto_cor_df_long,
      aes(x = Lag, y = Autocorrelation)
    ) +
      geom_line() +
      facet_wrap(~Parameter) +
      labs(
        title = plot_title
      )
  } else {
    # Plots!
    p <- ggplot(
      auto_cor_df_long,
      aes(x = Lag, y = Autocorrelation, colour = Parameter)
    ) +
      geom_line() +
      labs(
        title = plot_title
      )
  }
  p
}

# Save autocorrelation plots using ggplot2 style
saveAutoCorrelationPlots <- function(mcmc_data, save_dirs, plot_titles,
                                     n_seeds = length(mcmc_data),
                                     lag.max = NULL,
                                     gen_auto_corr_title = "Autocorrelation",
                                     plot_type = ".png",
                                     facet = FALSE,
                                     save_plots = FALSE) {

  # The filenames and plot titles for the plots
  auto_corr_file_names <- paste0(save_dirs, "/autocorrelation_plot", plot_type)
  auto_corr_titles <- paste0(plot_titles, ": ", gen_auto_corr_title)

  p_lst <- vector("list", n_seeds)

  for (i in 1:n_seeds) {


    # Plot title and filename
    curr_file_name <- auto_corr_file_names[i]
    curr_title <- auto_corr_titles[i]

    # Create the autocorrelation plot
    p_lst[[i]] <- autoCorrelationPlot(mcmc_data[[ii]], curr_title, lag.max = lag.max, facet = facet)

    if (save_plots) {
      ggsave(curr_file_name, plot = p_lst[[i]])
    }
  }
  # Return a list of plots
  p_lst
}

gewekePlot <- function(x,
                       frac_1 = 0.1,
                       frac_2 = 0.5,
                       n_bins = 20,
                       p_value_threshold = 0.05,
                       threshold_line_colour = "grey",
                       plt_title = "Geweke diagnostic plot") {

  # The preferred object type for interacting with coda functions
  x <- as.mcmc.list(x)

  # The vector of start iterations to calculate the Geweke statistic for
  start_iter_vec <- seq(
    from = start(x),
    to = (start(x) + end(x)) / 2,
    length = n_bins
  )

  # The matrix that will hold the Geweke stat
  geweke_mat <- matrix(nrow = length(start_iter_vec), ncol = nvar(x), dimnames = list(start_iter_vec, varnames(x)))

  for (n in 1:length(start_iter_vec)) {
    curr_geweke_diag <- geweke.diag(window(x, start = start_iter_vec[n]),
      frac1 = frac_1,
      frac2 = frac_2
    )

    geweke_mat[n, ] <- curr_geweke_diag[[1]]$z
  }

  # The 1.96 threshold for 0.05 significance on a standard normal distribution
  c_limit <- qnorm(1 - p_value_threshold / 2)

  # The variables to gather when moving from wide to long data (these are our
  # parameters)
  vars_to_gather <- varnames(x)

  # The data.frame we will plot (transform to long data to use the ggplot2
  # framework)
  geweke_df <- data.frame(Start_iteration = start_iter_vec) %>%
    cbind(geweke_mat) %>%
    tidyr::gather_("Parameter", "Geweke_statistic", vars_to_gather)

  # For this kind of plot I prefer unifrom axis, thus we find the y-axis limits
  y_limit <- max(c(c_limit, abs(geweke_df$Geweke_statistic)))

  # Plot the Geweke statistic for each parameter including a significance
  # threshold
  p <- ggplot(geweke_df, aes(x = Start_iteration, y = Geweke_statistic)) +
    geom_line() +
    geom_hline(yintercept = c_limit, linetype = "dashed", color = threshold_line_colour) +
    geom_hline(yintercept = -c_limit, linetype = "dashed", color = threshold_line_colour) +
    ylim(-y_limit, y_limit) +
    facet_wrap(~Parameter) +
    labs(
      x = "First iteration in segment",
      y = "Geweke's convergence dianostic",
      title = plt_title
    )

  p
}


saveGewekePlots <- function(mcmc_lst, save_dirs, plot_titles,
                            n_seeds = length(mcmc_lst),
                            gen_title = "Geweke plot",
                            plot_type = ".png",
                            save_plots = FALSE,
                            frac_1 = 0.1,
                            frac_2 = 0.5,
                            n_bins = 20,
                            p_value_threshold = 0.05,
                            threshold_line_colour = "grey") {


  # The filenames and plot titles for the plots
  file_names <- paste0(save_dirs, "/geweke_plot", plot_type)
  plot_titles <- paste0(plot_titles, ": ", gen_title)

  p_lst <- vector("list", n_seeds)

  for (i in 1:n_seeds) {
    x <- mcmc.list(as.mcmc(mcmc_lst[[ii]]))

    # Plot title and filename
    curr_file_name <- file_names[i]
    curr_title <- plot_titles[i]

    p_lst[[i]] <- gewekePlot(x,
      frac_1 = frac_1,
      frac_2 = frac_2,
      n_bins = n_bins,
      p_value_threshold = p_value_threshold,
      plt_title = curr_title,
      threshold_line_colour = threshold_line_colour
    )

    if (save_plots) {
      ggsave(curr_file_name, plot = p_lst[[i]])
    }
  }
  # Return a list of plots
  p_lst
}

gelmanValues <- function(x,
                         max_bins = 50,
                         confidence = 0.95,
                         transform = FALSE,
                         autoburnin = TRUE) {
  n_bins <- min(floor((niter(x) - 50) / thin(x)), max_bins)

  if (n_bins < 1) {
    stop("Insufficient iterations to produce Gelman-Rubin plot")
  }

  # Set the bin width within which to calculate the Shrinkage factor
  bin_width <- floor((niter(x) - 50) / n_bins)

  # The last iteration considered within each bin
  last_iter_vec <- c(
    seq(
      from = start(x) + 50 * thin(x),
      by = bin_width * thin(x),
      length = n_bins
    ),
    end(x)
  )


  confidence_threshold <- paste(50 * (confidence + 1), "%",
    sep = ""
  )

  # Array to hold the shrinkage factor (median and condifence threshold)
  shrink <- array(dim = c(n_bins + 1, nvar(x), 2))

  dimnames(shrink) <- list(
    last_iter_vec, varnames(x),
    c("median", confidence_threshold)
  )

  # Calculate the shrinkage factor for each bin of iterations
  for (i in 1:(n_bins + 1)) {
    shrink[i, , ] <- gelman.diag(window(x, end = last_iter_vec[i]),
      confidence = confidence,
      transform = transform,
      autoburnin = autoburnin,
      multivariate = FALSE
    )$psrf
  }

  all.na <- apply(
    is.na(shrink[, , 1, drop = FALSE]), 2,
    all
  )
  if (any(all.na)) {
    cat("Cannot compute Gelman & Rubin's diagnostic for any chain \n")
    cat("segments for variables", varnames(x)[all.na], "\n")
    cat("This indicates convergence failure\n")
    stop()
  }

  list("shrink" = shrink, "last_iter_vec" = last_iter_vec)
}


gelmanPlot <- function(x,
                       max_bins = 50,
                       confidence = 0.95,
                       transform = FALSE,
                       autoburnin = TRUE,
                       auto.layout = TRUE,
                       xlab = "last iteration in chain",
                       ylab = "shrink factor",
                       title = "Gelman-Rubin diagnostic plot") {
  x <- as.mcmc.list(x)

  gelman_obj <- gelmanValues(x,
    max_bins = max_bins,
    confidence = confidence,
    transform = transform,
    autoburnin = autoburnin
  )

  shrink <- gelman_obj$shrink
  last_iter_vec <- gelman_obj$last_iter_vec

  median_values <- shrink[, , "median"]
  threshold_values <- shrink[, , confidence_threhsold]

  vars_to_gather <- varnames(x)



  threshold_df <- data.frame(
    "Last_iter" = last_iter_vec,
    "Quantity" = confidence_threhsold
  ) %>%
    cbind(threshold_values) %>%
    tidyr::gather("Parameter", "Gelman_stat", vars_to_gather)


  median_df <- data.frame("Last_iter" = last_iter_vec, "Quantity" = "median") %>%
    cbind(median_values) %>%
    tidyr::gather("Parameter", "Gelman_stat", vars_to_gather)


  gelman_df <- rbind(median_df, threshold_df)

  p <- ggplot(gelman_df, aes(x = Last_iter, y = Gelman_stat, colour = Quantity, linetype = Quantity)) +
    geom_line() +
    scale_color_manual(values = c("black", "red")) +
    facet_wrap(~Parameter) +
    labs(
      x = xlab,
      y = ylab,
      title = title
    )

  p
}

input_arguments <- function() {
  option_list <- list(

    # Parent directory to MDI output directories
    optparse::make_option(c("-d", "--data_dir"),
      type = "character",
      default = ".",
      help = "Parent directory to MDI output directories [default= %default]",
      metavar = "character"
    ),

    # Generic beginning of MDI output files (I save them in the form out_seed_1.csv)
    optparse::make_option(c("--gen_input_file_name"),
      type = "character",
      default = "out_",
      help = "Generic beginning of MDI output files [default= %default]",
      metavar = "character"
    ),

    # File extension of MDI output files
    optparse::make_option(c("--file_ext"),
      type = "character",
      default = ".csv",
      help = "File extension of MDI output files [default= %default]",
      metavar = "character"
    ),

    # File type of plots to save
    optparse::make_option(c("--plot_type"),
      type = "character",
      default = ".png",
      help = "File type of plots to save [default= %default]",
      metavar = "character"
    ),

    # Instruction to save burn-in diagnostic plots
    optparse::make_option(c("--plot_burn_in"),
      type = "logical",
      default = TRUE,
      help = "Instruction to save burn-in diagnostic plots [default= %default]",
      metavar = "character"
    ),

    # Instruction to save autocorrelation diagnostic plots
    optparse::make_option(c("--plot_auto_corr"),
      type = "logical",
      default = TRUE,
      help = "Instruction to save autocorrelation diagnostic plots [default= %default]",
      metavar = "character"
    ),

    # Instruction to save Geweke diagnostic plots
    optparse::make_option(c("--plot_geweke"),
      type = "logical",
      default = TRUE,
      help = "Instruction to save Geweke diagnostic plots [default= %default]",
      metavar = "character"
    ),

    # Instruction to save Gelman-Rubin diagnostic plot
    optparse::make_option(c("--plot_gelman"),
      type = "logical",
      default = TRUE,
      help = "Instruction to save Gelman-Rubin diagnostic plot [default= %default]",
      metavar = "character"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}


# === Setup ====================================================================

args <- input_arguments()

# Set ggplot2 theme
theme_set(theme_bw())

# Example directory for output data
gen_mdi_dir <- args$data_dir #  "~/Documents/PhD/Year_1/Consensus_clustering/Yeast/Output_data/MDI/"

# List the sub directories
sub_dirs <- list.dirs(
  path = gen_mdi_dir,
  full.names = F,
  recursive = FALSE
)

# The generic file name before the seed and extension
gen_file_name <- args$gen_input_file_name # "out_"

# The file output type (this should never change)
file_ext <- args$file_ext # ".csv"

# One of .png or .pdf
plot_type <- args$plot_type # ".png"

# Instructions regarding plots to construct
plot_burn_in <- args$plot_burn_in # TRUE
plot_auto_corr <- args$plot_auto_corr # TRUE
plot_geweke <- args$plot_geweke # TRUE
plot_gelman <- args$plot_gelman # TRUE

# The file names will be in order of strings - 10 comes before 2. Remedy this.
numerical_order <- sub_dirs %>%
  extractNumericalFromString() %>%
  as.numeric() %>%
  order()

# Rearrange the sub directories to be intuitive
sub_dirs <- sub_dirs[numerical_order]
full_sub_dirs <- paste0(gen_mdi_dir, sub_dirs)

# This is the sub directory we will save our diagnostic plots to
new_sub_dirs <- "Convergence_diagnostics"
new_full_sub_dirs <- paste(full_sub_dirs, new_sub_dirs, sep = "/")

# Create these directories
# Set to an object to not print messages
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
  ess_p_lst <- saveESSPlots(my_data, new_full_sub_dirs, pretty_sub_dir_names,
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
  auto_p_lst <- saveAutoCorrelationPlots(new_data,
    new_full_sub_dirs,
    pretty_sub_dir_names,
    lag.max = NULL,
    gen_auto_corr_title = "Autocorrelation",
    plot_type = plot_type,
    facet = TRUE,
    save_plots = TRUE
  )
}

# auto_p_lst[[1]] / auto_p_lst[[10]] + plot_layout(guides = "collect")

# Plot the Geweke diagnostic statistic for each seed
if (plot_geweke) {
  geweke_p_lst <- saveGewekePlots(new_data,
    new_full_sub_dirs,
    pretty_sub_dir_names,
    save_plots = TRUE,
    threshold_line_colour = "red"
  )
}

# Plot the Gelman-Rubin shrinkage factor
gelmanPlot(new_data,
  xlab = "Last iteration in chain",
  ylab = "Shrinkage factor",
  title = "MDI: Gelman-Rubin diagnostic plot"
)
