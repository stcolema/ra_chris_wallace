
library(dplyr)
geweke_data <- read.csv("/Users/stephen/Desktop/Convergence/large_n_small_p_large_k_small_dm/simulation_1GewekeData.csv")


test_data <- geweke_data %>% 
  group_by(Chain) %>% 
  summarise(Shapiro_p_value = shapiro.test(Geweke_statistic)$p.value) %>% 
  mutate(Normal = Shapiro_p_value > 0.10)

match_ind <- match(geweke_data$Chain, test_data$Chain)
geweke_data$Converged <- test_data$Normal[match_ind]
