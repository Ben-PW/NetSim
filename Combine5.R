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

#### Importing data and functions ####
# This needs to be done before loading packages as some scripts load and detach packages
source('Error_functions_3.R')
source('Data_process.R')
source('Compute_NetStats.R')

#### Load environment ####

library(intergraph)
library(dplyr)

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

# If memory useage becomes an issue, skip creating the list, run the calculations step by step, then delete
# networks when no longer needed
# Or create multiple smaller lists, then delete networks when no longer needed and then run the other 
# list

error_levels <- c(0.1, 0.3, 0.5)

error_networks <- list(
  missingTies = randomMissingTies(datasets, error_levels),
  addedTies = randomAddedTies(datasets, error_levels),
  missingNodes = randomMissingNodes(datasets, error_levels),
  simpleAddedNodes = randomAddedNodesFixed(datasets, error_levels),
  DDAddedNodes = randomAddedNodesDD(datasets, error_levels)
)

##### Recalculate statistics on altered networks #####

error_networks <- map(error_networks, ~map(.x, ~map(.x, computeCentrality)))

metrics_error_networks <- map(error_networks, function(error_type_list) {
  imap(error_type_list, function(dataset_list, dataset_name) {
    map(dataset_list, function(glist) {
      computeMetrics(glist, paste0("E", dataset_name))
    })
  })
})


##### Compute bias #####

# This for if you want to calculate all at once, not recommended for large runs
centrality_bias_all <- imap(error_networks, function(error_type_list, error_type_name) { # error_type_name = missingTies, etc
  imap(error_type_list, function(dataset_list, dataset_name) { # dataset_name = Flo etc
    imap(dataset_list, function(perturbed_list, error_label) { # error_label = p03 etc
      computeNodeBias(
        original_sim = datasets[[dataset_name]],
        perturbed_sim = perturbed_list,
        name = paste0("E", dataset_name, "_", error_label),
        method = "pearson" # Remember to change if necessary
      )
    })
  })
})

##### Hashed out functions for applying error individually #####

# Missing ties
# this to calculate sections individually
#biasMissingTies <- purrr::imap(missingTies, function(error_levels, dataset_name) {
#  purrr::imap(error_levels, function(perturbed_list, error_label) {
#   computeNodeBias(
#      original_sim = datasets[[dataset_name]],
#      perturbed_sim = perturbed_list,
#      name = paste0("E", dataset_name, "_", error_label),
#      method = "pearson"
#    )
#  })
#})

# Added ties

#biasAddedTies <- purrr::imap(addedTies, function(error_levels, dataset_name) {
#  purrr::imap(error_levels, function(perturbed_list, error_label) {
#    computeNodeBias(
#      original_sim = datasets[[dataset_name]],
#      perturbed_sim = perturbed_list,
#      name = paste0("E", dataset_name, "_", error_label)
#    )
#  })
#})

# Missing nodes

#biasMissingNodes <- purrr::imap(missingNodes, function(error_levels, dataset_name) {
#  purrr::imap(error_levels, function(perturbed_list, error_label) {
#    computeNodeBias(
#      original_sim  = datasets[[dataset_name]],
#      perturbed_sim = perturbed_list,
#      name          = paste0("E", dataset_name, "_", error_label)
#    )
#  })
#})

# DD added nodes

#biasDDAddedNodes <- purrr::imap(DDAddedNodes, function(error_levels, dataset_name) {
#  purrr::imap(error_levels, function(perturbed_list, error_label) {
#    computeNodeBias(
#      original_sim  = datasets[[dataset_name]],
#      perturbed_sim = perturbed_list,
#      name          = paste0("E", dataset_name, "_", error_label)
#    )
#  })
#})


# Combine all results into a single dataframe
# metric_bias <- bind_rows(biasEFloMTies.01, biasEAc1MTies.01, biasEAc2MTies.01)

# Add in original full network ground truth metrics
# bias_GTadd <- left_join(metric_bias, GTAll, by = "id")

##### End of individual section #####

##### Extract data from error lists

centrality_bias_df <- imap_dfr(centrality_bias_all, function(error_type_list, error_type_name) {
  imap_dfr(error_type_list, function(dataset_list, dataset_name) {
    imap_dfr(dataset_list, function(df, error_label) {
      df %>%
        mutate(
          error_type  = error_type_name,
          dataset     = dataset_name,
          error_level = as.numeric(sub("p", ".", error_label))
        )
    })
  })
})

global_bias_df <- imap_dfr(metrics_error_networks, function(error_type_list, error_type_name) {
  imap_dfr(error_type_list, function(dataset_list, dataset_name) {
    imap_dfr(dataset_list, function(df, error_label) {
      df %>%
        mutate(
          error_type  = error_type_name,
          dataset     = dataset_name,
          error_level = as.numeric(sub("p", ".", error_label))
        )
    })
  })
})


########################################### Visualisations ######################################

library(tidyr)
library(ggplot2)

##### Prepare data for visualisation #####

# Global metrics
metrics_summary <- global_bias_df %>%
  pivot_longer(cols = c(density, dcent, clustering, APL),
               names_to = "Metric", values_to = "Value") %>%
  group_by(error_type, dataset, error_level, Metric) %>%
  summarise(
    mean_val = mean(Value, na.rm = TRUE),
    sd_val   = sd(Value, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(
    sd_lower = mean_val - sd_val,
    sd_upper = mean_val + sd_val
  )

# Node-level bias
bias_summary <- centrality_bias_df %>%
  pivot_longer(cols = starts_with(c("Degree", "Betweenness", "Closeness", "Eigenvector", "PageRank")),
               names_to = "Metric", values_to = "Correlation") %>%
  group_by(error_type, dataset, error_level, Metric) %>%
  summarise(
    mean_val = mean(Correlation, na.rm = TRUE),
    sd_val   = sd(Correlation, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(
    sd_lower = mean_val - 1.96 * sd_val,
    sd_upper = mean_val + 1.96 * sd_val
  )

##### Visualisation #####

plot_error_summary <- function(df, y_label, title_prefix, ylim = NULL) {
  ggplot(df, aes(x = error_level, y = mean_val,
                 color = dataset, group = dataset)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = sd_lower, ymax = sd_upper, fill = dataset),
                alpha = 0.2, colour = NA) +
    facet_grid(Metric ~ error_type, scales = "free_y") +
    labs(
      x = "Error level",
      y = y_label,
      title = paste(title_prefix, "by error type")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    ) +
    coord_cartesian(ylim = ylim)
}

# Global metric robustness
plot_error_summary(metrics_summary,
                   y_label = "Mean ± 1 SD",
                   title_prefix = "Global network metrics")

# Centrality robustness
plot_error_summary(bias_summary,
                   y_label = "Average correlation ± 1.96 SD",
                   title_prefix = "Centrality robustness",
                   ylim = c(-1, 1))
