####################################################################################################################

# Seeing if implementing these functions is easier in Igraph environment
# Turns out is is, by a country mile

#####################################################################################################################

library(igraph)

# Get a list of test IGraph objects
iFlo <- lapply(FloSim, intergraph::asIgraph)
iAc1 <- lapply(Ac1Sim, intergraph::asIgraph)
iAc2 <- lapply(Ac2Sim, intergraph::asIgraph)

iFlo[1]
iFlo[2]

# Coerce directed networks to symmetrical

undirect <- function(graph_list) {
  lapply(graph_list, function(g) {
    if (igraph::is.directed(g)) {
      g <- igraph::as.undirected(g, mode = "collapse")
    }
    g
  })
}

iAc1 <- undirect(iAc1)
iAc2 <- undirect(iAc2)

# Assign IDs for testing between node removal or addition situations

IDNodes <- function(graph_list){
  lapply(graph_list, function(g){
    V(g)$NodeID <- seq_len(vcount(g))
    g
  })
}

iFlo <- IDNodes(iFlo)
iAc1 <- IDNodes(iAc1)
iAc2 <- IDNodes(iAc2)
# V(iAc1[[1]])$NodeID

# Compute ground truth metrics

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

GTFlo <- computeMetrics(iFlo, "Flo")
GTAc1 <- computeMetrics(iAc1, "Ac1")
GTAc2 <- computeMetrics(iAc2, "Ac2")

GTAll <- bind_rows(GTFlo, GTAc1, GTAc2)

# Compute centrality scores and assign as attributes

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

iFlo <- computeCentrality(iFlo)
iAc1 <- computeCentrality(iAc1)
iAc2 <- computeCentrality(iAc2)

#V(iFlo[[1]])$Eigenvector
#V(iAc1[[100]])$Closeness

################################## Error functions for Igraph

# Randomly remove ties from the network
tieMissRand1 <- function(graph_list, missing_pct = 0.1) {
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

 EFlo1 <- tieMissRand1(iFlo, 0.2)

# Randomly remove nodes from the network
nodeMissRand1 <- function(graph_list, missing_pct = 0.1) {
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

# Eflo2 <- nodeMissRand1(iFlo, 0.2)
# V(Eflo2[[1]])$NodeID

# Function to randomly add ties
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

Eflo3 <- tieAddRand(iFlo, 0.2)

#edge_density(iFlo[[1]])
#edge_density(Eflo3[[1]])

# Add nodes with fixed number of random ties

addNodesRandFixed <- function(graph_list, error_pct = 0.1, edge_number = 2) {
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

Eflo4 <- addNodesRandFixed(iFlo)

iFlo[[1]]
Eflo4[[1]]
V(Eflo4[[1]])$NodeID

# Add nodes with ties based on original network degree distribution

# THIS FUNCTION IS NOT COMPLETE YET
nodeAddRandDD <- function(graph_list, error_pct = 0.1) {
  lapply(graph_list, function(g){
    
    n_old <- vcount(g)
    n_new <- round(error_pct * n_old)
    g <- add_vertices(g, n_new)
    new_vs <- which(is.na(V(g)$NodeID))
  })
    net <- network_list[[i]]  # Extract network object
    
    # Get number of original nodes
    num_original_nodes <- network.size(net)
    
    # Determine the number of new nodes to add
    num_new_nodes <- round(num_original_nodes * add_pct)
    
    if (num_new_nodes > 0) {
      # Add new nodes
      add.vertices(net, num_new_nodes)
      
      # Get total number of nodes after adding new ones
      num_total_nodes <- network.size(net)
      
      # Compute the indices of the newly added nodes
      new_node_indices <- (num_original_nodes + 1):num_total_nodes
      
      # Get degree distribution of the original nodes
      original_degrees <- degree(net)[1:num_original_nodes]  # Exclude new nodes from degree calc
      
      # Ensure there are existing nodes to connect to
      if (length(original_degrees) > 0) {
        for (new_node in new_node_indices) {
          # Determine the number of ties for the new node based on the degree distribution
          num_new_ties <- sample(original_degrees, 1)
          
          # Ensure we don't select more nodes than available
          num_new_ties <- min(num_new_ties, num_original_nodes)
          
          # Select random existing nodes to connect to
          selected_nodes <- sample(1:num_original_nodes, num_new_ties, replace = FALSE)
          
          # Add ties between the new node and selected existing nodes
          add.edges(net, tail = rep(new_node, length(selected_nodes)), head = selected_nodes)
        }
      }
    }
    
    # Store modified network
    network_list_added[[i]] <- net
  }
  
  return(network_list_added)
}

# Compute bias metrics
# This function is gross and should probably be changed somehow

compute_bias_metrics7 <- function(original_sim, perturbed_sim, name) {
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

EFlo1 <- computeCentrality(EFlo1)
FloSim1 <- compute_bias_metrics7(iFlo, EFlo1, "FloSim")

