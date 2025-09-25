#####################################################################################################
# Helper functions for simulation script
####################################################################################################

# Coerce directed networks to symmetrical

undirect <- function(graph_list) {
  lapply(graph_list, function(g) {
    if (igraph::is.directed(g)) {
      g <- igraph::as.undirected(g, mode = "collapse")
    }
    g
  })
}

# Assign IDs for testing between node removal or addition situations

IDNodes <- function(graph_list){
  lapply(graph_list, function(g){
    V(g)$NodeID <- seq_len(vcount(g))
    g
  })
}

# Return matching Node IDs

# This function did not return nodes in the same order, meaning correlations were 
# incorrect and artificially low

#matchNodes <- function(g1, g2) {
#  common_ids <- intersect(V(g1)$NodeID, V(g2)$NodeID)
#  list(
#    x = match(common_ids, V(g1)$NodeID),
#    y = match(common_ids, V(g2)$NodeID)
#  )
#}

matchNodes <- function(g1, g2) {
  common_ids <- intersect(V(g1)$NodeID, V(g2)$NodeID)
  
  # Order according to g1
  x <- match(common_ids, V(g1)$NodeID)
  
  # Get the same common_ids in that order then match in g2
  y <- match(V(g1)$NodeID[x], V(g2)$NodeID)
  
  list(x = x, y = y)
}

# Correlation wrapper for bias calculation

correlateNodes <- function(x, y, method = "pearson") {
  if (method == "pearson") {
    cor(x, y, use = "complete.obs")
  } else if (method == "spearman") {
    cor(rank(x), rank(y), method = "spearman", use = "complete.obs")
  } else {
    stop("Method must be 'pearson' or 'spearman'")
  }
}

##### Compute network level metrics #####

# I think I need to address the rounding here
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

##### Function to compute bias metrics when networks have been perturbed and may have missing
# or spuriously present nodes. Function will consider only nodes present in both ground truth AND
# error networks, using the NodeID vertex attribute as a key

#######################################################################################

# READ ME
# READ ME
# READ ME

# ADD A CALCULATION FOR KENDALL'S TAU RANK CORRELATION 
# AS SEEN IN
# https://link.springer.com/article/10.1140/epjb/e2007-00033-7

# Compute bias metrics

# This function is gross and should probably be changed somehow
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

centrality <- "Degree"

correlateDegree <- function(g1, g2, method){
  IDs <- matchNodes(original_sim, perturbed_sim)
  correlateNodes(
    V(g1)$Degree[IDs$x],
    V(g2)$Degree[IDs$y],
    method
  )
}

# I remember the problem here, you can't extract Igraph attributes using [[]], will need
# to find another method

#correlateCentrality <- function(g1, g2, centrality, method){
#  IDs <- matchNodes(g1, g2)
#  correlateNodes(
#    V(g1)[[centrality]][IDs$x],
#    V(g2)[[centrality]][IDs$y],
#    method
#  )
#}

correlateCentrality <- function(g1, g2, centrality, method){
  IDs <- matchNodes(g1, g2)
  x <- vertex_attr(g1, centrality, index = IDs$x)
  y <- vertex_attr(g2, centrality, index = IDs$y)
  cor(x, y, method = method, use = "complete.obs")
}

IDs <- matchNodes(iFlo[[100]], EFloMTies.05[[100]])
x <- vertex_attr(iFlo[[100]], "Degree", index = IDs$x)
y <- vertex_attr(EFloMTies.05[[100]], "Degree", index = IDs$y)
x
y

IDs <- matchNodes(iFlo[[100]], EFloMTies.05[[100]])

V(iFlo[[100]])$NodeID[IDs$x] == V(EFloMTies.05[[100]])$NodeID[IDs$y]

test1 <- iFlo[[100]]
test2 <- EFloMTies.05[[100]]

correlateCentrality(test1,test2, "Degree", "pearson")
correlateCentrality(test1,test2, "Closeness", "pearson")
correlateCentrality(test1,test2, "Eigenvector", "pearson")

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
