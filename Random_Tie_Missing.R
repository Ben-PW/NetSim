
##### Randomly remove ties #####
TieMissRand <- function(network_list, missing_pct = 0.1) {
  # Initialize an empty list to store networks with missing ties
  network_list_missing <- list()
  
  # Loop through each network in the list
  for (i in seq_along(network_list)) {
    
    # Extract the network object
    net <- as.matrix(network_list[[i]])

    # Get indices of unique ties (only upper triangle to avoid duplicate pairs)
    edge_indices <- which(upper.tri(net) & net != 0, arr.ind = TRUE)
    
    # Determine number of ties to remove
    num_edges <- nrow(edge_indices)
    num_remove <- ceiling(num_edges * missing_pct)  # Ensure at least one removal
    
    if (num_remove > 0 && num_edges > 0) {
      # Select edges to remove randomly
      edges_to_remove <- edge_indices[sample(1:num_edges, num_remove), , drop = FALSE]
      
      # Remove selected ties symmetrically
      for (j in 1:nrow(edges_to_remove)) {
        net[edges_to_remove[j, 1], edges_to_remove[j, 2]] <- 0
        net[edges_to_remove[j, 2], edges_to_remove[j, 1]] <- 0  # Maintain symmetry
      }
    }
    
    # Convert back to sna::network object and store the modified network
    network_list_missing[[i]] <- network::network(net, directed = FALSE)
  }
  
  return(network_list_missing)
}

##### Randomly remove nodes #####
NodeMissRand <- function(network_list, missing_pct = 0.1) {
  network_list_missing <- list()
  
  for (i in seq_along(network_list)) {
    net <- as.matrix(network_list[[i]])
    
    num_nodes <- nrow(net)
    num_remove <- ceiling(num_nodes * missing_pct)
    
    if (num_remove > 0 && num_nodes > 0) {
      nodes_to_remove <- sample(1:num_nodes, num_remove)
      
      net[nodes_to_remove, ] <- 0
      net[, nodes_to_remove] <- 0
    }
    
    network_list_missing[[i]] <- network::network(net, directed = FALSE)
  }
  
  return(network_list_missing)
}

##### Randomly add ties #####

TieAddRand <- function(network_list, add_pct = 0.1) {
  network_list_added <- list()
  
  for (i in seq_along(network_list)) {
    net <- as.matrix(network_list[[i]])
    
    num_nodes <- nrow(net)
    possible_edges <- which(upper.tri(net) & net == 0, arr.ind = TRUE)
    num_add <- ceiling(nrow(possible_edges) * add_pct)
    
    if (num_add > 0 && nrow(possible_edges) > 0) {
      edges_to_add <- possible_edges[sample(1:nrow(possible_edges), num_add), , drop = FALSE]
      
      for (j in 1:nrow(edges_to_add)) {
        net[edges_to_add[j, 1], edges_to_add[j, 2]] <- 1
        net[edges_to_add[j, 2], edges_to_add[j, 1]] <- 1  # Maintain symmetry
      }
    }
    
    network_list_added[[i]] <- network::network(net, directed = FALSE)
  }
  
  return(network_list_added)
}


##### Randomly add nodes #####
# This function will randomly add nodes to the network, with the number of ties added randomly
# based off the degrees of nodes in the original network

NodeAddRand_DD <- function(network_list, add_pct = 0.1) {
  network_list_added <- list()
  
  for (i in seq_along(network_list)) {
    net <- as.matrix(network_list[[i]])
    
    num_nodes <- nrow(net)
    num_add <- ceiling(num_nodes * add_pct)
    
    for (a in 1:num_add) {
      # Extract degree distribution
      degree_distribution <- network_list[[i]] %v% "Degree"
      
      # New node ID
      new_node_id <- nrow(net) + 1
      
      # Determine the number of ties for the new node
      num_new_ties <- sample(degree_distribution, 1)
      
      # Select nodes to connect to
      existing_nodes <- seq_len(nrow(net))
      nodes_to_connect <- sample(existing_nodes, num_new_ties, replace = FALSE)
      
      # Expand adjacency matrix
      net <- rbind(net, 0)
      net <- cbind(net, 0)
      
      # Add ties to the selected nodes
      for (node in nodes_to_connect) {
        net[new_node_id, node] <- 1
        net[node, new_node_id] <- 1
      }
    }
    
    # Store modified network
    network_list_added[[i]] <- network::network(net, directed = FALSE)
  }
  
  return(network_list_added)
}
