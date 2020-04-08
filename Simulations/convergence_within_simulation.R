#!/usr/env/bin/Rscript

# Various libraries
library(ggplot2)
library(mdiHelpR)
library(magrittr)
library(coda)
library(mcclust)
library(data.table)
library(patchwork)
library(stringr)

# Not used
# library(fitR) # devtools::install_github("sbfnk/fitR")

# My dewfault ggplot2 theme
setMyTheme(base_size = 14)

# Directory
analysis_dir <- "./PhD/Year_1/Consensus_inference/Consensus_inference_gen/Simulations/Simulation_results/Simulations/Single_dataset/base_case/simulation_1/Bayesian/"

# Read the files and sort by chain number (recogininsing as a numeric)
my_files <- list.files(analysis_dir, full.names = T) %>%
  stringr::str_sort(numeric = T)
file_order <- my_files %>% stringr::str_order(numeric = T)

n_files <- length(my_files)

samples <- list()
for (i in 1:n_files) {
  samples[[i]] <- fread(my_files[i], select = 1)
}


# Parameters regarding chain lengths and samples to keep
expected_length <- 1e6
burnin <- 1e5
thin_by <- 100

# Check all the chains have the same number of samples; stop analysis if not
fail <- sapply(samples, nrow) != expected_length
if (any(fail)) {
  print(paste("Chain(s)", paste(which(fail), collapse = " & "), "has less samples than others"))
  stop("Stopping analysis.")
}

# Apply burn in and thinning, and convert to mcmc type (for coda)
convergence_data <- list()
for (i in 1:n_files) {
  convergence_data[[i]] <- samples[[i]][seq(burnin, expected_length, thin_by), ] %>%
    # set_rownames(seq(1, expected_length, thin_by)) %>%
    set_colnames("Mass Parameter") %>%
    as.mcmc()
}

# Gelman-Rubin plot
p1 <- gelmanPlot(convergence_data) +
  labs(title = "Across chain convergence") +
  scale_x_continuous(label = function(x) {
    return(paste(x * thin_by + burnin))
  })

# Geweke plot for first chain
p2 <- gewekePlot(convergence_data[[1]],
  threshold_line_colour = "Red"
) +
  labs(
    title = "Within chain convergence",
    subtitle = "Chain 1"
  ) +
  scale_x_continuous(label = function(x) {
    return(paste(x * thin_by + burnin))
  }) # ,
# breaks = scales::pretty_breaks(3)            #### Note that 3 might be a bad choice
# )


# Get the geweke data for remaining chains
p_geweke <- lapply(
  convergence_data,
  gewekePlot
)

for (i in 2:n_files) {
  curr_df <- p_geweke[[i]]$data
  curr_df$Simulation <- paste0("Simulation ", i)

  if (i == 2) {
    g_df <- curr_df
  } else {
    g_df <- rbind(g_df, curr_df)
  }
}

g_df$Simulation <- g_df$Simulation %>%
  factor(levels = stringr::str_sort(unique(.), numeric = T))

# The 1.96 threshold for 0.05 significance on a standard normal distribution
c_limit <- stats::qnorm(1 - 0.05 / 2)

# Construct facet wrapped Geweke plots
y_limit <- max(c(c_limit, abs(g_df$Geweke_statistic)))
threshold_line_colour <- "red"
plt_title <- "Within chain convergence for remaining chains"
p_rest <- ggplot(g_df, aes(x = Start_iteration, y = Geweke_statistic)) +
  geom_line() +
  geom_hline(yintercept = c_limit, linetype = "dashed", color = threshold_line_colour) +
  geom_hline(yintercept = -c_limit, linetype = "dashed", color = threshold_line_colour) +
  ylim(-y_limit, y_limit) +
  facet_wrap(~Simulation) +
  labs(
    x = "First iteration in segment",
    y = "Geweke's convergence dianostic",
    title = plt_title,
    subtitle = "Mass parameter"
  ) +
  scale_x_continuous(
    label = function(x) {
      return(paste(x * thin_by + burnin))
    },
    breaks = scales::pretty_breaks(3)
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(3)
  )


# Use patchwork to create a grid of plots
(p1 + p2) / p_rest

ggplot_ratio <- 8.02 / 8.97
ggplot_width <- 12
ggsave("./PhD/Year_1/Consensus_inference/Consensus_inference_gen/Simulations/Simulation_pipe/example_convergence_within_sim.png",
  width = ggplot_width,
  heigh = ggplot_width * ggplot_ratio
)
