####################################################################################################################

# Error and helper functions currently stored here

#####################################################################################################################

# Dependenices 

library(igraph)
library(purrr)
library(dplyr)

################################## Error functions for Igraph

##### Random tie missingness #####

# Base function to randomly remove ties from the network

tieMissRand <- function(graph_list, error_pct = 0.1) {
  lapply(graph_list, function(g) {
    m <- ecount(g)
    k <- round(m * error_pct)
    if (m > 0 && k > 0) {
      # sample k random edges
      to_remove <- sample(E(g), k)
      # delete_edges returns a new graph
      g <- delete_edges(g, to_remove)
    }
    g
  })
}

# Wrapper to apply to network lists

randomMissingTies <- function(datasets, error_levels) {
  map(datasets, function(data) {
    map(error_levels, function(p) tieMissRand(data, error_pct = p)) %>%
      set_names(paste0("p", sub("\\.", "", error_levels)))
  })
}

##### Random node deletion #####

# Base function to remove nodes from the network

nodeMissRand <- function(graph_list, error_pct = 0.1) {
  lapply(graph_list, function(g) {
    n <- vcount(g)
    k <- round(n * error_pct)
    if (n > 0 && k > 0) {
      to_remove <- sample(V(g), k)
      # delete nodes and their edges
      g <- delete_vertices(g, to_remove)
    }
    g
  })
}

# Wrapper to apply to network lists

randomMissingNodes <- function(datasets, error_levels) {
  map(datasets, function(data) {
    map(error_levels, function(p) nodeMissRand(data, error_pct = p)) %>%
      set_names(paste0("p", sub("\\.", "", error_levels)))
  })
}

##### Random tie addition #####

# Base function to randomly add ties

tieAddRand <- function(graph_list, error_pct = 0.1) {
  lapply(graph_list, function(g) {
    m <- ecount(g)
    num_add <- round(m * error_pct)
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

# Wrapper to apply to network lists

randomAddedTies <- function(datasets, error_levels) {
  map(datasets, function(data) {
    map(error_levels, function(p) tieAddRand(data, error_pct = p)) %>%
      set_names(paste0("p", sub("\\.", "", error_levels)))
  })
}

##### Node addition functions #####

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

randomAddedNodesFixed <- function(datasets, error_levels, edge_number = 2) {
  map(datasets, function(data) {
    map(error_levels, function(p) nodeAddRand(data, error_pct = p, edge_number = edge_number)) %>%
      set_names(paste0("p", sub("\\.", "", error_levels)))
  })
}

# Add nodes with ties based on original network degree distribution

nodeAddRandDD <- function(graph_list, error_pct = 0.1) {
  lapply(graph_list, function(g) {

    n0    <- igraph::vcount(g)
    n_new <- round(n0 * error_pct)
    
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

randomAddedNodesDD <- function(datasets, error_levels) {
  map(datasets, function(data) {
    map(error_levels, function(p) nodeAddRandDD(data, error_pct = p)) %>%
      set_names(paste0("p", sub("\\.", "", error_levels)))
  })
}

# Cleanup

detach(package:dplyr)
detach(package:igraph)
detach(package:purrr)


