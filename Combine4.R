####################################################################################################

#Combine4.R is the updated version of Combine3.R. The aim of this script is to integrate
#a NodeID system, and update all the error and bias calculation functions to work with
#it, as well as update all error functions to rely on native sna functions. The idea is 
#that by matching nodes based on their ID, comparison of GT networks with error networks 
#that have had nodes added or removed is much more robust and simpler to implement 

#I succeeded in the above, but at what cost

####################################################################################################

detach(package:networkdata)

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

Ac1Sim <- undirect(iAc1)
Ac2Sim <- undirect(iAc2)

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
      # all‐pairs shortest paths; exclude Infs
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
      to_remove <- sample(E(g), k)
      # delete_edges returns a new graph
      g <- delete_edges(g, to_remove)
    }
    g
  })
}

# 10% random tie deletion
EFloMTies.01 <- tieMissRand(iFlo, 0.3)
EAc1MTies.01 <- tieMissRand(iAc1, 0.3)
EAc2MTies.01 <- tieMissRand(iAc2, 0.3)

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
        
        # Get their endpoint pairs and flatten for add_edges()
        ed_pairs <- ends(comp_g, to_add)
        edge_vec <- as.vector(t(ed_pairs))
        
        # Add 
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


##### Node level metrics #####

EFloMTies.01 <- computeCentrality(EFloMTies.01)
EAc1MTies.01 <- computeCentrality(EAc1MTies.01)
EAc2MTies.01 <- computeCentrality(EAc2MTies.01)

################################# Dataframe to compare centrality scores ##########################

##### Function to compute all relevant centrality metrics and store in dataframe #####

# This function is gross but the repeated intersect and match calls are actually really trivial computationally
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
biasEAc2MTies.01 <- computeNodeBias(iAc2, EAc1MTies.01, "Ac2")



# Combine all results into a single dataframe
metric_bias <- bind_rows(biasEFloMTies.01, biasEAc1MTies.01, biasEAc2MTies.01)

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

ggplot(bias_GTadd, aes(x = clustering, y = bet_robust_cor)) +
  geom_point(alpha = 0.6, color = "blue") +
  theme_minimal() 


