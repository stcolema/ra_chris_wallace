
library(magrittr)
library(ggplot2)
library(stringr)

mdiHelpR::setMyTheme()


data_dir <- "./Simulations/Results/"

scenarios <- list.dirs(data_dir, recursive = F, full.names = F)

scenarios_pretty <- scenarios %>% 
  str_replace_all("_", " ") %>% 
  str_to_sentence()

files <- paste0(data_dir, scenarios, "/all_results.csv")

n_scn <- length(scenarios)


my_data <- list()
for(i in 1:n_scn){
  my_data[[i]] <- read.csv(files[i], row.names = 1)
  my_data[[i]]$Scenario <- scenarios_pretty[i]
}

my_df <- do.call("rbind", my_data)

my_df$Scenario <- factor(my_df$Scenario, 
  levels = scenarios_pretty[c(8, 7, 1, 11, 5, 6, 9, 10, 2, 4, 3)])

model_labels <- my_df$Model %>% unique() %>% str_sort(numeric = T)
my_df$Model <- factor(my_df$Model, levels = model_labels[c(2:26, 1, 27)])

p1 <- my_df %>% 
  ggplot(aes(x = Model, y = ARI)) +
  geom_boxplot(colour = "black", fill = "#FDE725FF") +
  coord_flip() +
  facet_wrap(~Scenario)

ggsave(paste0(data_dir, "allmodel_performance_prediction.png"), width = 20, height = 16)


p2 <- my_df %>% 
  ggplot(aes(x = Model, y = Frobenius_norm)) +
  geom_boxplot(colour = "black", fill = "#FDE725FF") +
  coord_flip() +
  facet_wrap(~Scenario)

ggsave(paste0(data_dir, "allmodel_performance_uncertainty.png"), width = 20, height = 16)

library(dplyr)

na_check <- my_df %>% 
  group_by(Model, Scenario) %>% 
  summarise(NAs = sum(is.na(ARI)))
  # summarise(qs = function(x) 
  #   {if(sum(is.na(ARI) == nrow(ARI))){
  #     NA
  #     } 
  #     else {
  #       quantile(ARI, c(0.25, 0.75))
  #       }
  #     }, 
  #   prob = c(0.25, 0.75))

na_check[na_check$NAs == 99, ]

my_df %>% 
  dplyr::filter(Model != "Maximum likelihood (Mclust)" | Scenario != "No structure") %>% 
  group_by(Model, Scenario) %>% 
  summarise(Median = median(ARI))
# summarise(qs = quantile(ARI, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75))
