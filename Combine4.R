####################################################################################################

#Combine4.R is the updated version of Combine3.R. The aim of this script is to integrate
#a NodeID system, and update all the error and bias calculation functions to work with
#it, as well as update all error functions to rely on native sna functions. The idea is 
#that by matching nodes based on their ID, comparison of GT networks with error networks 
#that have had nodes added or removed is much more robust and simpler to implement 

#I succeeded in the above, but at what cost

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

####################################### Simulate networks #######################################

FloSim <- simulate(flo1, nsim = 1000)
Ac1Sim <- simulate(acct1, nsim = 1000)
Ac2Sim <- simulate(acct2, nsim = 1000)

rm(flo1, acct1, acct2)
################################ Ensure all networks are undirected ###########################

##### Function to ensure networks are undirected #####

undirect <- function(network_list) {
  lapply(network_list, function(net) {
    if (is.directed(net)) {
      # Convert to adjacency matrix, then back to an undirected network
      net <- as.network(as.matrix.network.adjacency(net), directed = FALSE)
    }
    
    return(net)
  })
}

##### Undirect networks #####

Ac1Sim <- undirect(Ac1Sim)
Ac2Sim <- undirect(Ac2Sim)

###################################### Assign ID to nodes ############################
##### Function to assign node ID #####
assign_node_ids <- function(network_list) {
  for (i in seq_along(network_list)) {
    network_list[[i]] %v% "NodeID" <- seq_len(network.size(network_list[[i]]))
  }
  return(network_list)
}

##### Assign IDs to GT nodes #####
FloSim <- assign_node_ids(FloSim)
Ac1Sim <- assign_node_ids(Ac1Sim)
Ac2Sim <- assign_node_ids(Ac2Sim)

############################### Compute Ground Truth metrics ####################################

##### Function to compute GT metrics for each simulation #####
compute_metrics <- function(network_list, name) {
  data.frame(
    id = paste0(name, seq_along(network_list)),
    density = sapply(network_list, network.density),
    dcent = centralization(network_list, degree),
    clustering = sapply(network_list, gtrans),
    size = sapply(network_list, network.size),
    avg_path_length = sapply(network_list, function(net) {
      gmean <- geodist(net)$gdist
      mean(gmean[is.finite(gmean)], na.rm = TRUE)
    })
  )
}

##### Compute GT metrics and combine #####

GTFlo <- compute_metrics(FloSim, "FloSim")
GTAc1 <- compute_metrics(Ac1Sim, "Ac1Sim")
GTAc2 <- compute_metrics(Ac2Sim, "Ac2Sim")

# Combine all into a single dataframe
GTAll <- bind_rows(GTFlo, GTAc1, GTAc2)

# Confirm statistic calculation and storage
cat("Ground truth network statistics calculated and stored in 'GTdf'.\n")

################################### Calculate node level metrics ################################

##### Function to compute and assign node level metrics #####
compute_centrality <- function(network_list) {
  convert_graph <- list()
  
  for(i in seq_along(network_list)) {
    convert_graph[[i]] <- intergraph::asIgraph(network_list[[i]])
    network_list[[i]] %v% "Degree" <- sna::degree(network_list[[i]], gmode = "graph")
    network_list[[i]] %v% "Betweenness" <- estimate_betweenness(convert_graph[[i]], cutoff = -1)
    network_list[[i]] %v% "Closeness" <- estimate_closeness(convert_graph[[i]], cutoff = -1)
    network_list[[i]] %v% "Eigenvector" <- eigen_centrality(convert_graph[[i]])$vector
    network_list[[i]] %v% "PageRank" <- page_rank(convert_graph[[i]])$vector
  }
  
  return(network_list)
}

##### Compute and assign node level metrics #####

library(igraph)

FloSim <- compute_centrality(FloSim)
Ac1Sim <- compute_centrality(Ac1Sim)
Ac2Sim <- compute_centrality(Ac2Sim)

detach(package:igraph)

######################################### Alter networks ########################################

source("Error_functions_SNA.R") 

##### Simple tie addition #####
# TieAddIterative is a simple and pretty efficient function, default use this because the 
# other ones scare me a bit

# CHECK BELOW IDEA TOMORROW #

# Would it be more memory efficient to do this with a network object which didn't have so many
# attribbutes? Like I could calculate initial centrality scores and store in a dataframe, 
# then delete everything except node id, add error, recalculate, store in another df, delete,
# rinse and repeat?

# 10% random tie addition
EFlo.01 <- TieAddRand6(FloSim, 0.1)
EAc1.01 <- TieAddIterative(Ac1Sim, 0.1)
EAc2.01 <- TieAddIterative(Ac2Sim, 0.1)
EAC2.01 <- TieAddRandNode(Ac2Sim, 0.1) # Takes quite long to run

# 20% random tie addition
EFlo.02 <- TieAddIterative(FloSim, 0.2)
EAc1.02 <- TieAddIterative(Ac1Sim, 0.2)
EAc2.02 <- TieAddIterative(Ac2Sim, 0.2)

##### Spurious nodes with 1 tie #####
EFlo.NA.01.01 <- NodeAddFixedEdges2(FloSim, 0.1, 1)
EAc1.NA.01.01 <- NodeAddFixedEdges2(Ac1Sim, 0.1, 1)
EAc2.NA.01.01 <- NodeAddFixedEdges2(Ac2Sim, 0.1, 1)

##### Spurious nodes with number of ties decided by network degree distribution



EAc1NodeAdd[[1]] %v% "NodeID"
#EFlo.02 <- TieAddIterative(FloSim, 0.1)
is.directed(Ac1Sim[[1]])
is.directed(FloSim[[1]])
EAc2.01[[1]]
Ac2Sim[[1]]
EFlo.01[[1]] %v% "NodeID"
EFlo.01[1]
FloSim[[1]]
EFlo.Test3[[1]]

intersect(Ac1Sim[[1]] %v% "NodeID", EAc1NodeAdd[[1]] %v% "NodeID")
##### TEST FOR NODE ADD FUNCTIONALITY #####

detach(package:igraph)
EfloTest <- NodeAddRand_DD(FloSim[1], 0.1)
FloTestdf <- compute_metrics(EFlo.03, "Test")

######################################### Store data ###########################################

##### Global network metrics #####
EFlo.01.df <- compute_metrics(EFlo.01, "EFlo")
EAc1.01.df <- compute_metrics(EAc1.01, "EAc1")
EAc2.01.df <- compute_metrics(EAc2.01, "EAc2")
EAc1NodeAdd.df <- compute_metrics(EAc1NodeAdd, "EAc1NodeAdd")

##### Node level metrics #####
library(igraph)

EFlo.01 <- compute_centrality(EFlo.01)
EAc1.01 <- compute_centrality(EAc1.01)
EAc2.01 <- compute_centrality(EAc2.01)
ENodeAdd <- compute_centrality(EAc1NodeAdd)

detach(package:igraph)


################################# Dataframe to compare centrality scores ##########################

# IMPORTANT: I have removed the code squaring the pearson correlation because I'm not 100% on 
# the justification for squaring it

##### Function to compute all relevant centrality metrics and store in dataframe #####

compute_bias_metrics7 <- function(original_sim, perturbed_sim, name) {
  data.frame(
    id = paste0(name, seq_along(original_sim)),
    
    deg_robust_cor = sapply(seq_along(original_sim), function(i) {
      # Match nodes by NodeID
      common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor((original_sim[[i]] %v% "Degree")[orig_indices], 
          (perturbed_sim[[i]] %v% "Degree")[pert_indices], 
          use = "complete.obs")
    }),
    
    bet_robust_cor = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor((original_sim[[i]] %v% "Betweenness")[orig_indices], 
          (perturbed_sim[[i]] %v% "Betweenness")[pert_indices], 
          use = "complete.obs")
    }),
    
    clo_robust_cor = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor((original_sim[[i]] %v% "Closeness")[orig_indices], 
          (perturbed_sim[[i]] %v% "Closeness")[pert_indices], 
          use = "complete.obs")
    }),
    
    eig_robust_cor = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor((original_sim[[i]] %v% "Eigenvector")[orig_indices], 
          (perturbed_sim[[i]] %v% "Eigenvector")[pert_indices], 
          use = "complete.obs")
    }),
    
    pg_robust_cor = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor((original_sim[[i]] %v% "PageRank")[orig_indices], 
          (perturbed_sim[[i]] %v% "PageRank")[pert_indices], 
          use = "complete.obs")
    }),
    
    deg_robust_srcc = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor(rank((original_sim[[i]] %v% "Degree")[orig_indices]), 
          rank((perturbed_sim[[i]] %v% "Degree")[pert_indices]), 
          method = "spearman", use = "complete.obs")
    }),
    
    bet_robust_srcc = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor(rank((original_sim[[i]] %v% "Betweenness")[orig_indices]), 
          rank((perturbed_sim[[i]] %v% "Betweenness")[pert_indices]), 
          method = "spearman", use = "complete.obs")
    }),
    
    clo_robust_srcc = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor(rank((original_sim[[i]] %v% "Closeness")[orig_indices]), 
          rank((perturbed_sim[[i]] %v% "Closeness")[pert_indices]), 
          method = "spearman", use = "complete.obs")
    }),
    
    eig_robust_srcc = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor(rank((original_sim[[i]] %v% "Eigenvector")[orig_indices]), 
          rank((perturbed_sim[[i]] %v% "Eigenvector")[pert_indices]), 
          method = "spearman", use = "complete.obs")
    }),
    
    pg_robust_srcc = sapply(seq_along(original_sim), function(i) {
      common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor(rank((original_sim[[i]] %v% "PageRank")[orig_indices]), 
          rank((perturbed_sim[[i]] %v% "PageRank")[pert_indices]), 
          method = "spearman", use = "complete.obs")
    })
  )
}

##### Compute bias #####

# Compute bias for all network simulations
bias_FloSim <- compute_bias_metrics7(FloSim, EFlo.01, "FloSim")
bias_FloSim1 <- compute_bias_metrics7(FloSim, EFlo.01, "FloSimTest")
bias_Ac1Sim <- compute_bias_metrics7(Ac1Sim, EAc1.01, "Ac1Sim")
bias_Ac2Sim <- compute_bias_metrics7(Ac2Sim, EAc2.01, "Ac2Sim")
bias_NodeAdd <- compute_bias_metrics7(Ac1Sim, ENodeAdd, "NodeTest")


# Combine all results into a single dataframe
metric_bias <- bind_rows(bias_FloSim, bias_Ac1Sim, bias_Ac2Sim)

# Add in original full network ground truth metrics
bias_GTadd <- left_join(metric_bias, GTAll, by = "id")


########################################### Visualisations ######################################

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

ggplot(bias_GTadd, aes(x = density, y = bet_robust_srcc)) +
  geom_point(alpha = 0.6, color = "blue") +
  theme_minimal() 


