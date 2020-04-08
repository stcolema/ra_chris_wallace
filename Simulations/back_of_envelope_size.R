
n_sim <- 100
n_scn <- 14
base_case_models_i <- 5.6
cyncical_correction <- 1.1
space_required <- (base_case_models_i * cyncical_correction) * n_sim * n_scn
space_required

# Probably much bigger due to large N small P scenarios

# if thin = 100
b1 <- 5.3
b2 <- 0.056
b1 - b2

saving <- (b1 - b2) * 100 * 14
saving

space_required - saving
