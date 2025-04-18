# Function to randomly remove ties #

library(igraph)
library(intergraph)

AddRandomTiesIgraph <- function(network_list, add_pct = 0.1) {
  network_list_added <- list()
  
  for (i in seq_along(network_list)) {
    sna_net <- network_list[[i]]  # Extract sna network
    
    # Convert to igraph
    graph <- intergraph::asIgraph(sna_net)
    
    # Get number of current edges
    num_existing_edges <- gsize(graph)
    
    # Determine the number of edges to add
    num_add <- ceiling(num_existing_edges * add_pct)
    
    # Get all possible node pairs
    num_nodes <- vcount(graph)
    all_possible_edges <- combn(V(graph), 2, simplify = FALSE)  # Unique node pairs
    
    # Filter only the pairs that are NOT connected
    unconnected_edges <- Filter(function(x) !are.connected(graph, x[1], x[2]), all_possible_edges)
    
    # Prevent adding more edges than available
    num_add <- min(num_add, length(unconnected_edges))
    
    if (num_add > 0) {
      # Randomly select new edges
      edges_to_add <- sample(unconnected_edges, num_add)
      
      # Convert selected edges to a vector format for add_edges()
      new_edges <- unlist(edges_to_add)
      
      # Add new edges to the graph
      graph <- add_edges(graph, new_edges)
    }
    
    # Convert back to sna
    network_list_added[[i]] <- intergraph::asNetwork(graph)
  }
  
  return(network_list_added)
}
