#####################################################################################################
# Helper functions for simulation script
####################################################################################################

##### Symmetrise and directed networks #####
undirect <- function(network_list) {
  lapply(network_list, function(net) {
    if (is.directed(net)) {
      # Convert to adjacency matrix, then back to an undirected network
      net <- as.network(as.matrix.network.adjacency(net), directed = FALSE)
    }
    
    return(net)
  })
}

##### Assign node IDs #####
assign_node_ids <- function(network_list) {
  for (i in seq_along(network_list)) {
    network_list[[i]] %v% "NodeID" <- seq_len(network.size(network_list[[i]]))
  }
  return(network_list)
}

##### Compute network level metrics #####
# Function to compute GT metrics for each simulation
compute_metrics <- function(network_list, name) {
  data.frame(
    id = paste0(name, seq_along(network_list)),
    edges = sapply(network_list, network.edgecount),
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

##### Function to compute all relevant centrality metrics and store in dataframe #####
compute_bias_metrics_old <- function(original_sim, perturbed_sim, name) {
  data.frame(
    id = paste0(name, seq_along(original_sim)),
    deg_robust_cor = sapply(seq_along(original_sim), function(i) {
      cor(original_sim[[i]] %v% "Degree", perturbed_sim[[i]] %v% "Degree", 
          use = "complete.obs")
    }),
    bet_robust_cor = sapply(seq_along(original_sim), function(i) {
      cor(original_sim[[i]] %v% "Betweenness", perturbed_sim[[i]] %v% "Betweenness", 
          use = "complete.obs")
    }),
    clo_robust_cor = sapply(seq_along(original_sim), function(i) {
      cor(original_sim[[i]] %v% "Closeness", perturbed_sim[[i]] %v% "Closeness", 
          use = "complete.obs")
    }),
    eig_robust_cor = sapply(seq_along(original_sim), function(i) {
      cor(original_sim[[i]] %v% "Eigenvector", perturbed_sim[[i]] %v% "Eigenvector", 
          use = "complete.obs")
    }),
    pg_robust_cor = sapply(seq_along(original_sim), function(i) {
      cor(original_sim[[i]] %v% "PageRank", perturbed_sim[[i]] %v% "PageRank", 
          use = "complete.obs")
    }),
    deg_robust_srcc = sapply(seq_along(original_sim), function(i) {
      cor(rank(original_sim[[i]] %v% "Degree"), rank(perturbed_sim[[i]] %v% "Degree"), 
          method = "spearman", use = "complete.obs")
    }),
    bet_robust_srcc = sapply(seq_along(original_sim), function(i) {
      cor(rank(original_sim[[i]] %v% "Betweenness"), rank(perturbed_sim[[i]] %v% "Betweenness"), 
          method = "spearman", use = "complete.obs")
    }),
    clo_robust_srcc = sapply(seq_along(original_sim), function(i) {
      cor(rank(original_sim[[i]] %v% "Closeness"), rank(perturbed_sim[[i]] %v% "Closeness"), 
          method = "spearman", use = "complete.obs")
    }),
    eig_robust_srcc = sapply(seq_along(original_sim), function(i) {
      cor(rank(original_sim[[i]] %v% "Eigenvector"), rank(perturbed_sim[[i]] %v% "Eigenvector"), 
          method = "spearman", use = "complete.obs")
    }),
    pg_robust_srcc = sapply(seq_along(original_sim), function(i) {
      cor(rank(original_sim[[i]] %v% "PageRank"), rank(perturbed_sim[[i]] %v% "PageRank"), 
          method = "spearman", use = "complete.obs")
    })
  )
}

###################################### Let's fucking try this again ###########################

##### Function to compute bias metrics when networks have been perturbed and may have missing
# or spuriously present nodes. Function will consider only nodes present in both ground truth AND
# error networks, using the NodeID vertex attribute as a key

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

#common_nodes <- seq_len(min(network.size(FloSim[[1]]), network.size(EFlo.03[[1]])))
#cor(rank(get.vertex.attribute(FloSim[[1]], "Degree"))[common_nodes], 
#      rank(get.vertex.attribute(EFlo.03[[1]], "Degree"))[common_nodes], 
#      method = "spearman", use = "complete.obs")

# Check we get the same results
#bias_FloSim <- compute_bias_metrics_old(FloSim, EFlo.01, "FloSim")
#bias_FloSim1 <- compute_bias_metrics7(FloSim, EFlo.01, "FloSim")
#identical(bias_FloSim, bias_FloSim1) # Came back TRUE!!!
