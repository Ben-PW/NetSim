####################################################################################################################

# Error and helper functions currently stored here

#####################################################################################################################

################################## Error functions for Igraph

# Randomly remove ties from the network

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



# Randomly remove nodes from the network

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

# Add nodes with fixed number of random ties

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


# Add nodes with ties based on original network degree distribution

NodeAddRandDD_igraph <- function(graph_list, add_pct = 0.1) {
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




