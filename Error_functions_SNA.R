################################# Remove ties randomly ##################################

TieMissRand2 <- function(network_list, missing_pct = 0.1) {

  network_list_missing <- list()
  

  for (i in seq_along(network_list)) {
  
    net <- network_list[[i]]
    
    # Get all existing edges
    edge_list <- as.matrix.network.edgelist(net)
    
    # Determine the number of ties to remove
    num_edges <- nrow(edge_list)
    num_remove <- round(num_edges * missing_pct)
    
    if (num_remove > 0 && num_edges > 0) {
      # Select random edges to remove
      edges_to_remove <- edge_list[sample(1:num_edges, num_remove), , drop = FALSE]
      
      # Remove selected edges
      for (j in 1:nrow(edges_to_remove)) {
        delete.edges(net, get.edgeIDs(net, v = edges_to_remove[j, 1], alter = edges_to_remove[j, 2]))
      }
    }
    
    # Store modified network
    network_list_missing[[i]] <- net
  }
  
  return(network_list_missing)
}

###################################### Remove nodes randomly ###################################

NodeMissRand3 <- function(network_list, missing_pct = 0.1) {
  network_list_missing <- list()
  
  for (i in seq_along(network_list)) {
    net <- network_list[[i]]  # Extract network object
    
    num_nodes <- network.size(net)  # Get number of nodes
    num_remove <- round(num_nodes * missing_pct)  # Calculate nodes to remove
    
    if (num_remove > 0 && num_nodes > 0) {
      # Get current node names
      node_names <- network.vertex.names(net)
      
      # Select node names to remove randomly
      nodes_to_remove <- sample(node_names, num_remove, replace = FALSE)
      
      # Convert node names to numeric indices
      node_indices <- which(node_names %in% nodes_to_remove)
      
      # Remove selected nodes using indices
      delete.vertices(net, node_indices)
    }
    
    # Store modified network
    network_list_missing[[i]] <- net
  }
  
  return(network_list_missing)
}

######################################## Add ties randomly #####################################

#### Other functions #####

TieAddRand4 <- function(network_list, add_pct = 0.1) {
  network_list_added <- list()
  
  for (i in seq_along(network_list)) {
    net <- network_list[[i]]  # Extract network object
    
    # Get all current edges
    existing_edges <- as.matrix.network.edgelist(net)
    num_existing_edges <- nrow(existing_edges)  # Number of edges in the original network
    
    # Get all possible node pairs (excluding existing edges)
    all_nodes <- network.vertex.names(net)  # Node names
    node_indices <- seq_along(all_nodes)    # Corresponding numeric indices
    possible_edges <- combn(node_indices, 2, simplify = FALSE)  # Generate all possible pairs
    
    # Remove already existing edges
    existing_edges_set <- apply(existing_edges, 1,function(x) paste(sort(x),collapse = "-"))
    possible_edges <- possible_edges[!sapply(possible_edges,function(x) paste(sort(x),collapse = "-")) 
                                     %in% existing_edges_set]
    
    # Determine the number of ties to add based on **current edges, not possible edges**
    num_add <- round(num_existing_edges * add_pct)
    
    # Prevent adding more edges than available possible edges
    num_add <- min(num_add, length(possible_edges))
    
    if (num_add > 0 && length(possible_edges) > 0) {
      # Select random pairs to add edges
      edges_to_add <- sample(possible_edges, num_add)
      
      # Convert to numeric indices for add.edges()
      tail_nodes <- sapply(edges_to_add, `[[`, 1)
      head_nodes <- sapply(edges_to_add, `[[`, 2)
      
      # Add new edges
      add.edges(net, tail = tail_nodes, head = head_nodes)
    }
    
    # Store modified network
    network_list_added[[i]] <- net
  }
  
  return(network_list_added)
}

########################### A more efficient tie adding function #######################
# This thing is disgusting but might be th eonly option for really large networks

TieAddRand5 <- function(network_list, add_pct = 0.1) {
  network_list_added <- list()
  
  for (i in seq_along(network_list)) {
    net <- network_list[[i]]  # Extract network object
    
    # Get existing edges (numeric indices)
    existing_edges <- as.matrix.network.edgelist(net)
    num_existing_edges <- nrow(existing_edges)  # Number of edges in the original network
    
    # Create a set of existing edges for fast lookup
    existing_edge_set <- setNames(rep(TRUE, num_existing_edges), apply(existing_edges, 1, function(x) paste(sort(x), collapse = "-")))
    
    # Determine the number of ties to add based on current edges
    num_add <- ceiling(num_existing_edges * add_pct)
    
    # Get node indices
    num_nodes <- network.size(net)
    
    if (num_add > 0 && num_nodes > 1) {
      # Efficiently sample **missing** edges instead of generating all pairs
      edges_to_add <- list()
      attempts <- 0
      
      while (length(edges_to_add) < num_add && attempts < num_add * 10) {
        # Randomly select two distinct nodes
        new_edge <- sample(1:num_nodes, 2, replace = FALSE)
        edge_key <- paste(sort(new_edge), collapse = "-") # makes search more efficient but 
                                                # also unsuitable for directed networks
        # Check if this edge already exists
        if (!edge_key %in% existing_edge_set) {
          edges_to_add[[length(edges_to_add) + 1]] <- new_edge
          existing_edge_set[edge_key] <- TRUE  # Add to existing edge set to prevent duplication
        }
        
        attempts <- attempts + 1
      }
      
      if (length(edges_to_add) > 0) {
        # Convert list to vector format for add.edges()
        tail_nodes <- sapply(edges_to_add, `[[`, 1)
        head_nodes <- sapply(edges_to_add, `[[`, 2)
        
        # Add new edges
        add.edges(net, tail = tail_nodes, head = head_nodes)
      }
    }
    
    # Store modified network
    network_list_added[[i]] <- net
  }
  
  return(network_list_added)
}

# More robust version of TieAddRand5
TieAddRand6 <- function(network_list, add_pct = 0.1, max_attempts = 10000) {
  network_list_added <- list()
  
  for (i in seq_along(network_list)) {
    net <- network_list[[i]]  # Extract network object
    
    # Get existing edges (numeric indices)
    existing_edges <- as.matrix.network.edgelist(net)
    num_existing_edges <- nrow(existing_edges)  # Number of edges in the original network
    
    # Create a set of existing edges for fast lookup
    existing_edge_set <- setNames(rep(TRUE, num_existing_edges), apply(existing_edges, 1, function(x) paste(sort(x), collapse = "-")))
    
    # Determine the number of ties to add based on current edges
    num_add <- round(num_existing_edges * add_pct)
    
    # Get node indices
    num_nodes <- network.size(net)
    
    if (num_add > 0 && num_nodes > 1) {
      # Efficiently sample **missing** edges instead of generating all pairs
      edges_to_add <- list()
      attempts <- 0
      
      while (length(edges_to_add) < num_add && attempts < num_add * max_attempts) {
        # Randomly select two distinct nodes
        new_edge <- sample(1:num_nodes, 2, replace = FALSE)
        edge_key <- paste(sort(new_edge), collapse = "-") # makes search more efficient but 
                                                        # also unsuitable for directed networks
        # Check if this edge already exists
        if (!edge_key %in% existing_edge_set) {
          edges_to_add[[length(edges_to_add) + 1]] <- new_edge
          existing_edge_set[edge_key] <- TRUE  # Add to existing edge set to prevent duplication
        }
        
        attempts <- attempts + 1
      }
      
      # Warn if not all requested edges were added
      if (length(edges_to_add) < num_add) {
        warning(paste(
          "Only added", length(edges_to_add), "of", num_add,
          "requested ties in network", i,
          "- possibly too dense. Increase `max_attempts` if needed."
        ))
      }
      
      if (length(edges_to_add) > 0) {
        # Convert list to vector format for add.edges()
        tail_nodes <- sapply(edges_to_add, `[[`, 1)
        head_nodes <- sapply(edges_to_add, `[[`, 2)
        
        # Add new edges
        add.edges(net, tail = tail_nodes, head = head_nodes)
      }
    }
    
    # Store modified network
    network_list_added[[i]] <- net
  }
  
  return(network_list_added)
}
########################### More understandable tie adding function #############################
# This will not do well for larger, denser networks

TieAddIterative <- function(network_list, add_pct = 0.1) {
  network_list_added <- list()
  
  for (i in seq_along(network_list)) {
    net <- network_list[[i]]  # Extract network object
    
    # Get list of nodes and total existing ties
    num_nodes <- network.size(net)
    existing_edges <- as.matrix.network.edgelist(net)
    num_existing_edges <- nrow(existing_edges)
    
    # Determine the number of ties to add
    num_add <- round(num_existing_edges * add_pct)
    
    # Initialize counter for added edges
    ties_added <- 0
    
    # Loop until we've added the required number of ties
    while (ties_added < num_add) {
      # Randomly pick a node
      node <- sample(1:num_nodes, 1)
      
      # Get indices of nodes it is NOT connected to
      current_ties <- get.neighborhood(net, node, type = "both")
      possible_targets <- setdiff(1:num_nodes, c(node, current_ties))
      
      # If there are available nodes to connect to, add a tie
      if (length(possible_targets) > 0) {
        target_node <- sample(possible_targets, 1)
        add.edge(net, tail = node, head = target_node)
        ties_added <- ties_added + 1
      }
    }
    
    # Store modified network
    network_list_added[[i]] <- net
  }
  
  return(network_list_added)
}


############################## Add nodes to network randomly with x ties ################

NodeAddFixedEdges2 <- function(network_list, error_fraction = 0.1, edge_number = 2) {
  network_list_added <- list()
  
  for (i in seq_along(network_list)) {
    net <- network_list[[i]]  # Extract network object
    
    # Get number of nodes before adding new ones
    num_original_nodes <- network.size(net)
    
    # Determine the number of nodes to add
    num_new_nodes <- round(num_original_nodes * error_fraction)
    
    if (num_new_nodes > 0) {
      # Add new nodes (without relying on names)
      add.vertices(net, num_new_nodes)
      
      # Get total number of nodes after adding new ones
      num_total_nodes <- network.size(net)
      
      # Compute the indices of the newly added nodes
      new_node_indices <- (num_original_nodes + 1):num_total_nodes
      
      # Get numeric indices of existing nodes (excluding new ones)
      existing_node_indices <- 1:num_original_nodes  # These are the original nodes
      
      # Ensure there are existing nodes to connect to
      if (length(existing_node_indices) > 0) {
        for (new_node in new_node_indices) {
          # Select random existing nodes for connections
          selected_nodes <- sample(existing_node_indices, min(edge_number, length(existing_node_indices)), replace = FALSE)
          
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

NodeAddRandDD <- function(network_list, add_pct = 0.1) {
  network_list_added <- list()
  
  for (i in seq_along(network_list)) {
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
