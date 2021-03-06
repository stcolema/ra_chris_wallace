#!/usr/bin/env Rscript

# Function to save histogram of phi parameter from MDI over MCMC iterations

plot_phi_histograms <- function(phis, file_path, start_index, eff_n_iter) {
  loc_dir <- paste0(file_path, "Phi_histograms/")
  dir.create(loc_dir, showWarnings = FALSE)

  for (i in 1:length(phis)) {
    for (j in 1:ncol(phis[[i]])) {
      curr_phi <- colnames(phis[[i]])[j]
      plot_name <- paste0(loc_dir, curr_phi, "_plot", plot_type)

      # Find which datasets are used here
      dataset_indices <- readr::parse_number(curr_phi)

      # The following steps do not work if we have more than 9 datasets (which is not relevant to me)
      # This extracts the first dataset index (i.e. from 67 takes 6)
      dataset_index_1 <- dataset_indices[1] %>%
        "/"(10) %>%
        floor(.)

      # This extracts the secodn dataset index
      dataset_index_2 <- dataset_indices[1] %>% "%%"(10)

      dataset_1 <- dataset_names[dataset_index_1]
      dataset_2 <- dataset_names[dataset_index_2]

      histogram_plot_name <- paste0(loc_dir, curr_phi, "_histogram_plot", plot_type)

      histogram_title <- paste(curr_phi, ": histogram plot")

      histogram_subtitle <- paste(
        "Iterations",
        (burn + thin),
        "through",
        n_iter,
        "for",
        dataset_1,
        "and",
        dataset_2
      )

      ggplot(data = phis[[1]][start_index:eff_n_iter, ], aes_string(x = curr_phi)) +
        geom_histogram() +
        labs(
          title = histogram_title,
          subtitle = histogram_subtitle
        )

      ggsave(histogram_plot_name)
    }
  }
}
