####################################################################################################

#Combine4.R is the updated version of Combine3.R. The aim of this script is to integrate
#a NodeID system, and update all the error and bias calculation functions to work with
#it, as well as update all error functions to rely on native sna functions. The idea is 
#that by matching nodes based on their ID, comparison of GT networks with error networks 
#that have had nodes added or removed is much more robust and simpler to implement 

#I have also refactored the code to rely on the igraph package. This has made all functions more robust 
#and far simpler

####################################################################################################



################################# Randomisation test ###########################################
# Alter this seed to check simulated networks are actually changing
set.seed(2345)

########################################### Begin ##############################################
#### Load environment ####
library(intergraph)
library(dplyr)


#### Importing data ####

source('Data_process.R')
source('Compute_NetStats.R')

detach(package:networkdata)

####################################### Simulate networks #######################################

iFlo <- simulate(flo1, nsim = 1000)
iAc1 <- simulate(acct1, nsim = 1000)
iAc2 <- simulate(acct2, nsim = 1000)

rm(flo1, acct1, acct2)

detach(package:statnet)
detach(package:sna)

library(igraph)
library(intergraph)

iFlo <- lapply(iFlo, intergraph::asIgraph)
iAc1 <- lapply(iAc1, intergraph::asIgraph)
iAc2 <- lapply(iAc2, intergraph::asIgraph)
################################ Ensure all networks are undirected ###########################

##### Function to ensure networks are undirected #####

undirect <- function(graph_list) {
  lapply(graph_list, function(g) {
    if (igraph::is.directed(g)) {
      g <- igraph::as.undirected(g, mode = "collapse")
    }
    g
  })
}

##### Undirect networks #####

iAc1 <- undirect(iAc1)
iAc2 <- undirect(iAc2)

###################################### Assign ID to nodes ############################
##### Function to assign node ID #####
IDNodes <- function(graph_list){
  lapply(graph_list, function(g){
    V(g)$NodeID <- seq_len(vcount(g))
    g
  })
}

##### Assign IDs to GT nodes #####
iFlo <- IDNodes(iFlo)
iAc1 <- IDNodes(iAc1)
iAc2 <- IDNodes(iAc2)

############################### Compute Ground Truth metrics ####################################

##### Function to compute GT metrics for each simulation #####
computeMetrics <- function(graph_list, name) {
  data.frame(
    id = paste0(name, seq_along(graph_list)),
    
    density = sapply(graph_list, function(g) {
      edge_density(g, loops = FALSE)
    }),
    
    dcent = sapply(graph_list, function(g) {
      cent <- centralization.degree(g, mode = "all", normalized = TRUE)
      cent$centralization
    }),
    
    clustering = sapply(graph_list, function(g) {
      transitivity(g, type = "global")
    }),
    
    size = sapply(graph_list, vcount),
    
    APL = sapply(graph_list, function(g) {
      # all pairs shortest paths - exclude inf
      dist_mat <- distances(g, mode = "all")
      mean(dist_mat[is.finite(dist_mat)], na.rm = TRUE)
    })
  )
}

##### Compute GT metrics and combine #####

GTFlo <- computeMetrics(iFlo, "Flo")
GTAc1 <- computeMetrics(iAc1, "Ac1")
GTAc2 <- computeMetrics(iAc2, "Ac2")

# Combine all into a single dataframe
GTAll <- bind_rows(GTFlo, GTAc1, GTAc2)



################################### Calculate node level metrics ################################

##### Function to compute and assign node level metrics #####
computeCentrality <- function(graph_list){
  lapply(graph_list, function(g){
    V(g)$Degree <- degree(g, mode = "all")
    V(g)$Betweenness <- estimate_betweenness(g, cutoff = -1)
    V(g)$Closeness <- estimate_closeness(g, cutoff = -1)
    V(g)$Eigenvector <- eigen_centrality(g)$vector
    V(g)$PageRank <- page_rank(g)$vector
    g
  })
}

##### Compute and assign node level metrics #####

iFlo <- computeCentrality(iFlo)
iAc1 <- computeCentrality(iAc1)
iAc2 <- computeCentrality(iAc2)

######################################### Alter networks ########################################

source("Error_functions_3.R") 

##### Simple tie removal #####

# Tie removal function
tieMissRand <- function(graph_list, missing_pct = 0.1) {
  lapply(graph_list, function(g) {
    m <- ecount(g)
    k <- round(m * missing_pct)
    if (m > 0 && k > 0) {
      # sample k random edges
      to_remove <- sample(E(g), k, replace = F)
      
      g <- delete_edges(g, to_remove)
    }
    g
  })
}

# 10% random tie deletion
EFloMTies.01 <- tieMissRand(iFlo, 0.1)
EAc1MTies.01 <- tieMissRand(iAc1, 0.1)
EAc2MTies.01 <- tieMissRand(iAc2, 0.1)

EFloMTies.02 <- tieMissRand(iFlo, 0.2)
EAc1MTies.02 <- tieMissRand(iAc1, 0.2)
EAc2MTies.02 <- tieMissRand(iAc2, 0.2)

EFloMTies.03 <- tieMissRand(iFlo, 0.3)
EAc1MTies.03 <- tieMissRand(iAc1, 0.3)
EAc2MTies.03 <- tieMissRand(iAc2, 0.3)

EFloMTies.04 <- tieMissRand(iFlo, 0.4)
EAc1MTies.04 <- tieMissRand(iAc1, 0.4)
EAc2MTies.04 <- tieMissRand(iAc2, 0.4)

EFloMTies.05 <- tieMissRand(iFlo, 0.5)
EAc1MTies.05 <- tieMissRand(iAc1, 0.5)
EAc2MTies.05 <- tieMissRand(iAc2, 0.5)

##### Simple tie addition #####

# Tie addition function
tieAddRand <- function(graph_list, add_pct = 0.1) {
  lapply(graph_list, function(g) {
    m <- ecount(g)
    num_add <- round(m * add_pct)
    if (m > 0 && num_add > 0) {
      # Build complement graph 
      comp_g <- complementer(g, loops = FALSE)
      
      missing_edges <- ecount(comp_g)
      num_add       <- min(num_add, missing_edges) #Don't add more edges than is possible
      if (num_add > 0) {
        
        to_add <- sample(E(comp_g), num_add)
        
        # Get endpoints pairs and flatten for add_edges()
        ed_pairs <- ends(comp_g, to_add)
        edge_vec <- as.vector(t(ed_pairs))
        
        g <- add_edges(g, edge_vec)
      }
    }
    g
  })
}

# 10% random tie addition
EFloATies.01 <- tieAddRand(iFlo, 0.1)
EAc1ATies.01 <- tieAddRand(iAc1, 0.1)
EAc2ATies.01 <- tieAddRand(iAc2, 0.1)

##### Simple node deletion #####

# Node removal function
nodeMissRand <- function(graph_list, missing_pct = 0.1) {
  lapply(graph_list, function(g) {
    n <- vcount(g)
    k <- round(n * missing_pct)
    if (n > 0 && k > 0) {
      to_remove <- sample(V(g), k)
      # delete nodes and their edges
      g <- delete_vertices(g, to_remove)
    }
    g
  })
}

# 10% random node deletion
EFloMNodes.01 <- nodeMissRand(iFlo, 0.1)
EAc1MNodes.01 <- nodeMissRand(iAc1, 0.1)
EAc2MNodes.01 <- nodeMissRand(iAc2, 0.1)

##### Simple node addition #####

# Simple node addition function
nodeAddRand <- function(graph_list, error_pct = 0.1, edge_number = 2) {
  lapply(graph_list, function(g) {
    
    n_old <- vcount(g)
    n_new <- round(n_old * error_pct)
    if (n_new <= 0) return(g) 
    g <- add_vertices(g, n_new)
    new_vs <- which(is.na(V(g)$NodeID))
    
    # for each new v, sample 'edge_number' distinct old nodes
    edge_vec <- integer(0)
    for (v in new_vs) {
      k <- min(edge_number, n_old)
      targets <- sample(seq_len(n_old), k, replace = FALSE)
      # rbind(tail, head) then flatten
      pairs <- rbind(rep(v, k), targets)
      edge_vec <- c(edge_vec, as.vector(pairs))
    }
    
    g <- add_edges(g, edge_vec)
    g
  })
}

# 10% random node addition
EFloANodes.01 <- nodeAddRand(iFlo, 0.1, 2)
EAc1ANodes.01 <- nodeAddRand(iAc1, 0.1, 2)
EAc2ANodes.01 <- nodeAddRand(iAc2, 0.1, 2)

##### Degree distribution node addition #####

# Function to add nodes with N ties determined by sampling degree distribution
NodeAddRandDD <- function(graph_list, add_pct = 0.1) {
  lapply(graph_list, function(g) {
    
    n0    <- igraph::vcount(g)
    n_new <- round(n0 * add_pct)
    
    if (n_new > 0) {
      g      <- igraph::add_vertices(g, n_new)
      new_vs <- seq(n0 + 1, igraph::vcount(g))
      deg0   <- igraph::degree(g, v = seq_len(n0), mode = "all") # vector of node degrees
      
      for (v in new_vs) {
        k <- sample(deg0, 1) # sample existing node degrees
        k <- min(k, n0)
        if (k > 0) {
          peers <- sample(seq_len(n0), k)
          g     <- igraph::add_edges(g, c(rbind(v, peers))) # randomly add ties
          # option here to weight sampling by degree, preferential attachment style
        }
      }
    }
    g
  })
}

# 10% random node addition with degree distribution sampling
EFloANodesDD.01 <- nodeAddRand(iFlo, 0.1)
EAc1ANodesDD.01 <- nodeAddRand(iAc1, 0.1)
EAc2ANodesDD.01 <- nodeAddRand(iAc2, 0.1)

######################################### Store data ###########################################

# Get the network metrics of the perturbed networks

##### Global network metrics #####
EFloMTies.01.df <- computeMetrics(EFloMTies.01, "EFlo")
EAc1MTies.01.df <- computeMetrics(EAc1MTies.01, "EAc1")
EAc2MTies.01.df <- computeMetrics(EAc2MTies.01, "EAc2")

EFloMTies.02.df <- computeMetrics(EFloMTies.02, "EFlo")
EAc1MTies.02.df <- computeMetrics(EAc1MTies.02, "EAc1")
EAc2MTies.02.df <- computeMetrics(EAc2MTies.02, "EAc2")

EFloMTies.03.df <- computeMetrics(EFloMTies.03, "EFlo")
EAc1MTies.03.df <- computeMetrics(EAc1MTies.03, "EAc1")
EAc2MTies.03.df <- computeMetrics(EAc2MTies.03, "EAc2")

EFloMTies.04.df <- computeMetrics(EFloMTies.04, "EFlo")
EAc1MTies.04.df <- computeMetrics(EAc1MTies.04, "EAc1")
EAc2MTies.04.df <- computeMetrics(EAc2MTies.04, "EAc2")

EFloMTies.05.df <- computeMetrics(EFloMTies.05, "EFlo")
EAc1MTies.05.df <- computeMetrics(EAc1MTies.05, "EAc1")
EAc2MTies.05.df <- computeMetrics(EAc2MTies.05, "EAc2")


##### Node level metrics #####

EFloMTies.01 <- computeCentrality(EFloMTies.01)
EAc1MTies.01 <- computeCentrality(EAc1MTies.01)
EAc2MTies.01 <- computeCentrality(EAc2MTies.01)

EFloMTies.02 <- computeCentrality(EFloMTies.02)
EAc1MTies.02 <- computeCentrality(EAc1MTies.02)
EAc2MTies.02 <- computeCentrality(EAc2MTies.02)

EFloMTies.03 <- computeCentrality(EFloMTies.03)
EAc1MTies.03 <- computeCentrality(EAc1MTies.03)
EAc2MTies.03 <- computeCentrality(EAc2MTies.03)

EFloMTies.04 <- computeCentrality(EFloMTies.04)
EAc1MTies.04 <- computeCentrality(EAc1MTies.04)
EAc2MTies.04 <- computeCentrality(EAc2MTies.04)

EFloMTies.05 <- computeCentrality(EFloMTies.05)
EAc1MTies.05 <- computeCentrality(EAc1MTies.05)
EAc2MTies.05 <- computeCentrality(EAc2MTies.05)


################################# Dataframe to compare centrality scores ##########################

##### Function to compute all relevant centrality metrics and store in dataframe #####

# This function is gross but the repeated intersect and match calls are trivial computationally
# what really needs to change is to separate these out into separate functions for each metric and have them
# return a value which I can store in a dataframe normally, not by creating one

# IMPORTANT: I have removed the code squaring the pearson correlation because I'm not 100% on 
# the justification for squaring it

computeNodeBias <- function(original_sim, perturbed_sim, name) {
  data.frame(
    id = paste0(name, seq_along(original_sim)),
    
    deg_robust_cor = sapply(seq_along(original_sim), function(i) {
      # Match nodes by NodeID
      common_ids <- intersect(V(original_sim[[i]])$NodeID, V(perturbed_sim[[i]])$NodeID)
      orig_indices <- match(common_ids, V(original_sim[[i]])$NodeID)
      pert_indices <- match(common_ids, V(perturbed_sim[[i]])$NodeID)
      
      cor((V(original_sim[[i]])$Degree)[orig_indices], 
          (V(perturbed_sim[[i]])$Degree)[pert_indices], 
          use = "complete.obs")
    }),
    
    bet_robust_cor = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(V(original_sim[[i]])$NodeID, V(perturbed_sim[[i]])$NodeID)
      orig_indices <- match(common_ids, V(original_sim[[i]])$NodeID)
      pert_indices <- match(common_ids, V(perturbed_sim[[i]])$NodeID)
      
      cor((V(original_sim[[i]])$Betweenness)[orig_indices], 
          (V(perturbed_sim[[i]])$Betweenness)[pert_indices], 
          use = "complete.obs")
    }),
    
    clo_robust_cor = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(V(original_sim[[i]])$NodeID, V(perturbed_sim[[i]])$NodeID)
      orig_indices <- match(common_ids, V(original_sim[[i]])$NodeID)
      pert_indices <- match(common_ids, V(perturbed_sim[[i]])$NodeID)
      
      cor((V(original_sim[[i]])$Closeness)[orig_indices], 
          (V(perturbed_sim[[i]])$Closeness)[pert_indices], 
          use = "complete.obs")
    }),
    
    eig_robust_cor = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(V(original_sim[[i]])$NodeID, V(perturbed_sim[[i]])$NodeID)
      orig_indices <- match(common_ids, V(original_sim[[i]])$NodeID)
      pert_indices <- match(common_ids, V(perturbed_sim[[i]])$NodeID)
      
      cor((V(original_sim[[i]])$Eigenvector)[orig_indices], 
          (V(perturbed_sim[[i]])$Eigenvector)[pert_indices], 
          use = "complete.obs")
    }),
    
    pg_robust_cor = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(V(original_sim[[i]])$NodeID, V(perturbed_sim[[i]])$NodeID)
      orig_indices <- match(common_ids, V(original_sim[[i]])$NodeID)
      pert_indices <- match(common_ids, V(perturbed_sim[[i]])$NodeID)
      
      cor((V(original_sim[[i]])$PageRank)[orig_indices], 
          (V(perturbed_sim[[i]])$PageRank)[pert_indices], 
          use = "complete.obs")
    }),
    
    deg_robust_srcc = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(V(original_sim[[i]])$NodeID, V(perturbed_sim[[i]])$NodeID)
      orig_indices <- match(common_ids, V(original_sim[[i]])$NodeID)
      pert_indices <- match(common_ids, V(perturbed_sim[[i]])$NodeID)
      
      cor(rank((V(original_sim[[i]])$Degree)[orig_indices]), 
          rank((V(perturbed_sim[[i]])$Degree)[pert_indices]), 
          method = "spearman", use = "complete.obs")
    }),
    
    bet_robust_srcc = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(V(original_sim[[i]])$NodeID, V(perturbed_sim[[i]])$NodeID)
      orig_indices <- match(common_ids, V(original_sim[[i]])$NodeID)
      pert_indices <- match(common_ids, V(perturbed_sim[[i]])$NodeID)
      
      cor(rank((V(original_sim[[i]])$Betweenness)[orig_indices]), 
          rank((V(perturbed_sim[[i]])$Betweenness)[pert_indices]), 
          method = "spearman", use = "complete.obs")
    }),
    
    clo_robust_srcc = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(V(original_sim[[i]])$NodeID, V(perturbed_sim[[i]])$NodeID)
      orig_indices <- match(common_ids, V(original_sim[[i]])$NodeID)
      pert_indices <- match(common_ids, V(perturbed_sim[[i]])$NodeID)
      
      cor(rank((V(original_sim[[i]])$Closeness)[orig_indices]), 
          rank((V(perturbed_sim[[i]])$Closeness)[pert_indices]), 
          method = "spearman", use = "complete.obs")
    }),
    
    eig_robust_srcc = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(V(original_sim[[i]])$NodeID, V(perturbed_sim[[i]])$NodeID)
      orig_indices <- match(common_ids, V(original_sim[[i]])$NodeID)
      pert_indices <- match(common_ids, V(perturbed_sim[[i]])$NodeID)
      
      cor(rank((V(original_sim[[i]])$Eigenvector)[orig_indices]), 
          rank((V(perturbed_sim[[i]])$Eigenvector)[pert_indices]), 
          method = "spearman", use = "complete.obs")
    }),
    
    pg_robust_srcc = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(V(original_sim[[i]])$NodeID, V(perturbed_sim[[i]])$NodeID)
      orig_indices <- match(common_ids, V(original_sim[[i]])$NodeID)
      pert_indices <- match(common_ids, V(perturbed_sim[[i]])$NodeID)
      
      cor(rank((V(original_sim[[i]])$PageRank)[orig_indices]), 
          rank((V(perturbed_sim[[i]])$PageRank)[pert_indices]), 
          method = "spearman", use = "complete.obs")
    })
  )
}

##### Compute bias #####

# Compute bias for random tie missingness simulations
biasEFloMTies.01 <- computeNodeBias(iFlo, EFloMTies.01, "Flo")
biasEAc1MTies.01 <- computeNodeBias(iAc1, EAc1MTies.01, "Ac1")
biasEAc2MTies.01 <- computeNodeBias(iAc2, EAc2MTies.01, "Ac2")

biasEFloMTies.02 <- computeNodeBias(iFlo, EFloMTies.02, "Flo")
biasEAc1MTies.02 <- computeNodeBias(iAc1, EAc1MTies.02, "Ac1")
biasEAc2MTies.02 <- computeNodeBias(iAc2, EAc2MTies.02, "Ac2")

biasEFloMTies.03 <- computeNodeBias(iFlo, EFloMTies.03, "Flo")
biasEAc1MTies.03 <- computeNodeBias(iAc1, EAc1MTies.03, "Ac1")
biasEAc2MTies.03 <- computeNodeBias(iAc2, EAc2MTies.03, "Ac2")

biasEFloMTies.04 <- computeNodeBias(iFlo, EFloMTies.04, "Flo")
biasEAc1MTies.04 <- computeNodeBias(iAc1, EAc1MTies.04, "Ac1")
biasEAc2MTies.04 <- computeNodeBias(iAc2, EAc2MTies.04, "Ac2")

biasEFloMTies.05 <- computeNodeBias(iFlo, EFloMTies.05, "Flo")
biasEAc1MTies.05 <- computeNodeBias(iAc1, EAc1MTies.05, "Ac1")
biasEAc2MTies.05 <- computeNodeBias(iAc2, EAc2MTies.05, "Ac2")



# Combine all results into a single dataframe
metric_bias <- bind_rows(biasEFloMTies.01, biasEAc1MTies.01, biasEAc2MTies.01)

# Add in original full network ground truth metrics
bias_GTadd <- left_join(metric_bias, GTAll, by = "id")


########################################### Visualisations ######################################

# Placeholder visualisations for now 

library(ggplot2)
library(reshape2)

bias_long <- melt(metric_bias, 
                  id.vars = "id", 
                  variable.name = "Centrality", 
                  value.name = "Robustness")

ggplot(data = bias_long,
       aes(x = Centrality,
           y = Robustness)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(bias_GTadd, aes(x = clustering, y = bet_robust_cor)) +
  geom_point(alpha = 0.6, color = "blue") +
  theme_minimal() 

########################################### Bigger visualisation ##########################

#Combine metric dataframes

library(tidyr)

# Add error level column to each df
EFloMTies.01.df$error <- 0.1
EAc1MTies.01.df$error <- 0.1
EAc2MTies.01.df$error <- 0.1

EFloMTies.02.df$error <- 0.2
EAc1MTies.02.df$error <- 0.2
EAc2MTies.02.df$error <- 0.2

EFloMTies.03.df$error <- 0.3
EAc1MTies.03.df$error <- 0.3
EAc2MTies.03.df$error <- 0.3

EFloMTies.04.df$error <- 0.4
EAc1MTies.04.df$error <- 0.4
EAc2MTies.04.df$error <- 0.4

EFloMTies.05.df$error <- 0.5
EAc1MTies.05.df$error <- 0.5
EAc2MTies.05.df$error <- 0.5

#######################################################

metrics_all <- bind_rows(
  EFloMTies.01.df, EFloMTies.02.df, EFloMTies.03.df, EFloMTies.04.df, EFloMTies.05.df,
  EAc1MTies.01.df, EAc1MTies.02.df, EAc1MTies.03.df, EAc1MTies.04.df, EAc1MTies.05.df,
  EAc2MTies.01.df, EAc2MTies.02.df, EAc2MTies.03.df, EAc2MTies.04.df, EAc2MTies.05.df
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
biasEFloMTies.01$error <- 0.1; biasEAc1MTies.01$error <- 0.1; biasEAc2MTies.01$error <- 0.1
biasEFloMTies.02$error <- 0.2; biasEAc1MTies.02$error <- 0.2; biasEAc2MTies.02$error <- 0.2
biasEFloMTies.03$error <- 0.3; biasEAc1MTies.03$error <- 0.3; biasEAc2MTies.03$error <- 0.3
biasEFloMTies.04$error <- 0.4; biasEAc1MTies.04$error <- 0.4; biasEAc2MTies.04$error <- 0.4
biasEFloMTies.05$error <- 0.5; biasEAc1MTies.05$error <- 0.5; biasEAc2MTies.05$error <- 0.5


bias_all <- bind_rows(
  biasEFloMTies.01, biasEFloMTies.02, biasEFloMTies.03, biasEFloMTies.04, biasEFloMTies.05,
  biasEAc1MTies.01, biasEAc1MTies.02, biasEAc1MTies.03, biasEAc1MTies.04, biasEAc1MTies.05,
  biasEAc2MTies.01, biasEAc2MTies.02, biasEAc2MTies.03, biasEAc2MTies.04, biasEAc2MTies.05
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


