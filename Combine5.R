##################################################################################################

# Combine 5 is the updated version of combine 4. The aim of this version is to eliminate a lot
# of spaghetti code through the use of nested lists and some function updates, making the code far
# more robust and scalable

# This may not be fun

##################################################################################################

################################# Randomisation test ###########################################
# Alter this seed to check simulated networks are actually changing
set.seed(2345)

########################################### Begin ##############################################

#### Load environment ####

library(intergraph)
library(dplyr)

#### Importing data and functions ####
source('Error_functions_3.R')
source('Data_process.R')
source('Compute_NetStats.R')

####################################### Simulate networks #######################################

iFlo <- simulate(flo1, nsim = 1000)
iAc1 <- simulate(acct1, nsim = 1000)
iAc2 <- simulate(acct2, nsim = 1000)

# Inserting a larger network

ibook <- replicate(1000, network.copy(books), simplify = FALSE)

rm(flo1, acct1, acct2)

detach(package:statnet)
detach(package:sna)

library(igraph)
library(intergraph)

##### Define datasets #####

datasets <- list(
  Flo = iFlo,
  Ac1 = iAc1,
  Ac2 = iAc2,
  ibook = ibook
)

##### Convert all to igraph objects #####

datasets <- lapply(datasets, function(netlist) {
  lapply(netlist, intergraph::asIgraph)
})

##### Ensure all are undirected #####

datasets <- lapply(datasets, undirect)

##### Assign ID to nodes #####

datasets <- lapply(datasets, IDNodes)

##### Compute Ground Truth metrics #####

GTAll <- Map(computeMetrics, datasets, names(datasets)) %>%
  bind_rows()

##### Calculate node level metrics #####

datasets <- lapply(datasets, computeCentrality)

##### Alter networks #####

library(purrr)

# Simple random tie removal 
missingTies <- randomMissingTies(datasets, error_levels)

# Simple random tie addition
addedTies <- randomAddedTies(datasets, error_levels)

# Simple random node deletion 
missingNodes <- randomMissingNodes(datasets, error_levels)

# Simple random node addition (fixed edge number)
simpleAddedNodes <- randomAddedNodesFixed(datasets, error_levels)

# Random degree distribution node addition
DDAddedNodes <- randomAddedNodesDD(datasets, error_levels)

##### Store data #####

# Global network metrics

metricsMissingTies <- purrr::imap(missingTies, function(error_level_list, dataset_name) {
  purrr::map(error_level_list, function(glist) {
    computeMetrics(glist, paste0("E", dataset_name))
  })
})

# Node level metrics 

missingTies <- purrr::map(missingTies, function(error_level_list) {
  purrr::map(error_level_list, computeCentrality)
})

##### Compute bias #####

biasMissingTies <- purrr::imap(missingTies, function(error_levels, dataset_name) {
  purrr::imap(error_levels, function(perturbed_list, error_label) {
    computeNodeBias(
      original_sim  = datasets[[dataset_name]],
      perturbed_sim = perturbed_list,
      name          = paste0("E", dataset_name, "_", error_label)
    )
  })
})


# Combine all results into a single dataframe
metric_bias <- bind_rows(biasEFloMTies.01, biasEAc1MTies.01, biasEAc2MTies.01)

# Add in original full network ground truth metrics
bias_GTadd <- left_join(metric_bias, GTAll, by = "id")


########################################### Visualisations ######################################

# Placeholder visualisations for now 

library(ggplot2)
library(reshape2)
library(tidyr)

# Add error level column to each df
EFloMTies.01.df$error <- 0.1
EAc1MTies.01.df$error <- 0.1
EAc2MTies.01.df$error <- 0.1
EibookMTies.01.df$error <- 0.1

EFloMTies.02.df$error <- 0.2
EAc1MTies.02.df$error <- 0.2
EAc2MTies.02.df$error <- 0.2
EibookMTies.02.df$error <- 0.2

EFloMTies.03.df$error <- 0.3
EAc1MTies.03.df$error <- 0.3
EAc2MTies.03.df$error <- 0.3
EibookMTies.03.df$error <- 0.3

EFloMTies.04.df$error <- 0.4
EAc1MTies.04.df$error <- 0.4
EAc2MTies.04.df$error <- 0.4
EibookMTies.04.df$error <- 0.4

EFloMTies.05.df$error <- 0.5
EAc1MTies.05.df$error <- 0.5
EAc2MTies.05.df$error <- 0.5
EibookMTies.05.df$error <- 0.5

#######################################################

metrics_all <- bind_rows(
  EFloMTies.01.df, EFloMTies.02.df, EFloMTies.03.df, EFloMTies.04.df, EFloMTies.05.df,
  EAc1MTies.01.df, EAc1MTies.02.df, EAc1MTies.03.df, EAc1MTies.04.df, EAc1MTies.05.df,
  EAc2MTies.01.df, EAc2MTies.02.df, EAc2MTies.03.df, EAc2MTies.04.df, EAc2MTies.05.df,
  EibookMTies.01.df, EibookMTies.02.df, EibookMTies.03.df, EibookMTies.04.df, EibookMTies.05.df
)


metrics_all <- metrics_all %>% select(-size)


metrics_long <- metrics_all %>%
  tidyr::pivot_longer(cols = c(density, dcent, clustering, APL),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(Network = substr(id, 1, 4))   # "EFlo", "EAc1", "EAc2"

# Summarise mean and SD by metric, error, and network
metrics_summary <- metrics_long %>%
  group_by(Metric, error, Network) %>%
  summarise(
    mean_val = mean(Value, na.rm = TRUE),
    sd_val   = sd(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sd_lower = mean_val - sd_val,
    sd_upper = mean_val + sd_val
  )

######################################## Plot global metrics ########################################

ggplot(metrics_summary, aes(x = error, y = mean_val, color = Network, group = Network)) +
  geom_line() +
  geom_ribbon(aes(ymin = sd_lower, ymax = sd_upper, fill = Network), alpha = 0.2, colour = NA) +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(x = "Missingness proportion",
       y = "Average value ± 1 SD",
       title = "Global network metrics under increasing tie missingness") +
  theme_minimal()

################################# Same for node level

# Combine centrality robustness 

# Add error level 
biasEFloMTies.01$error <- 0.1; biasEAc1MTies.01$error <- 0.1; biasEAc2MTies.01$error <- 0.1; biasEibookMTies.01$error <- 0.1
biasEFloMTies.02$error <- 0.2; biasEAc1MTies.02$error <- 0.2; biasEAc2MTies.02$error <- 0.2; biasEibookMTies.02$error <- 0.2
biasEFloMTies.03$error <- 0.3; biasEAc1MTies.03$error <- 0.3; biasEAc2MTies.03$error <- 0.3; biasEibookMTies.03$error <- 0.3
biasEFloMTies.04$error <- 0.4; biasEAc1MTies.04$error <- 0.4; biasEAc2MTies.04$error <- 0.4; biasEibookMTies.04$error <- 0.4
biasEFloMTies.05$error <- 0.5; biasEAc1MTies.05$error <- 0.5; biasEAc2MTies.05$error <- 0.5; biasEibookMTies.05$error <- 0.5


bias_all <- bind_rows(
  biasEFloMTies.01, biasEFloMTies.02, biasEFloMTies.03, biasEFloMTies.04, biasEFloMTies.05,
  biasEAc1MTies.01, biasEAc1MTies.02, biasEAc1MTies.03, biasEAc1MTies.04, biasEAc1MTies.05,
  biasEAc2MTies.01, biasEAc2MTies.02, biasEAc2MTies.03, biasEAc2MTies.04, biasEAc2MTies.05,
  biasEibookMTies.01, biasEibookMTies.02, biasEibookMTies.03, biasEibookMTies.04, biasEibookMTies.05
)


bias_long <- bias_all %>%
  tidyr::pivot_longer(cols = -c(id, error),
                      names_to = "Metric", values_to = "Correlation") %>%
  mutate(Network = substr(id, 1, 3))

# Summarise mean ± 1 SD by metric, error, and network
bias_summary <- bias_long %>%
  group_by(Metric, error, Network) %>%
  summarise(
    mean_val = mean(Correlation, na.rm = TRUE),
    sd_val   = sd(Correlation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sd_lower = mean_val - 1.96*sd_val,
    sd_upper = mean_val + 1.96*sd_val
  )

######################################## Plot centrality robustness ########################################

ggplot(bias_summary, aes(x = error, y = mean_val, color = Network, group = Network)) +
  geom_line() +
  geom_ribbon(aes(ymin = sd_lower, ymax = sd_upper, fill = Network), alpha = 0.2, colour = NA) +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(x = "Missingness proportion",
       y = "Average correlation ± 1.96 SD",
       title = "Centrality robustness under increasing tie missingness") +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal()

