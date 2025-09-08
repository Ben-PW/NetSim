#################################################################################################

# Function to simulate the spotlight effect on a network

# This function will take a network as an argument. Eventually the seed node will be specifiable
# as an argument passed into the function. Edge observation probabilities will be a function of
# each node's geodesic distance from the seed node.

# Differing decay parameters will be possible, such as exponential vs linear decay 

# Differing combination probabilities will be possible for incident nodes to determine tie
# Observation probability, such as mean vs harmonic mean

#################################################################################################

##### Version 1 #####

# This version is just to annotate the geodesic distance of each other node in the network from
# the seed

library(igraph)

# Spotlight seed inverse-distance 
spotlight <- function(g) {
  stopifnot(is.igraph(g), !is_directed(g))
  
  # Random seed node
  seed <- sample(V(g), 1)
  
  # Geodesic distances
  d <- as.numeric(distances(g, v = seed, to = V(g), mode = "all"))
  
  # Inverse distance: 1/d, with 0 for self and unreachable nodes
  inv_d <- ifelse(is.infinite(d), NA,
                  ifelse(d == 0, 1, 1/d)) # set seed as 1 so average between seed and
                                          # neighbours = 1, meaning ties maximally observed
  
  # Add as node attribute
  V(g)$dist_seed <- inv_d
  
  # Store which seed was used
  graph_attr(g, "seed") <- V(g)$NodeID[seed]
  
  return(g)
}

iAc2[[1]]

test <- spotlight(Ac2Sim[[1]]) # function actually works

V(test)$dist_seed # distances are being accurately stored

diameter(test) # inverse geodesic reflects actual path lengths

plot(test)

##### Version 2 #####

# This version is to implement the tie scoring system based on node distance from seed

# Spotlight annotator: node- and edge-level distances
spotlight <- function(g, 
                      seed = NULL,
                      combiner = c("mean", "min", "max", "sum", "harmonic")) {
  
  stopifnot(is.igraph(g), !is_directed(g))
  combiner <- match.arg(combiner)
  
  # Pick a seed if not supplied
  if (is.null(seed)) {
    seed <- sample(V(g), 1)     
  } else {
    seed <- which(V(g)$NodeID == seed)
  }
  
  # Compute geo distances
  d <- distances(g, v = seed, to = V(g), mode = "all")
  d <- as.numeric(d)
  
  # Inverse geo distance: seed = 1, finite d>0 = 1/d, unreachable = NA
  inv_d <- ifelse(is.infinite(d), NA,
                  ifelse(d == 0, 1, 1/d))
  
  V(g)$dist_seed <- inv_d
  
  # Edge-level distances via combiner
  ends_mat <- ends(g, E(g), names = FALSE)
  
  di <- V(g)$dist_seed
  
  # Figure out which of these is most suitable later, for now it was pretty easy to 
  # just include them all 
  
  d_edge <- switch(combiner,
                   mean = rowMeans(cbind(di[ends_mat[,1]], di[ends_mat[,2]]), na.rm = TRUE), # mean of node distances
                   min = pmin(di[ends_mat[,1]], di[ends_mat[,2]], na.rm = TRUE), # min node distance
                   max = pmax(di[ends_mat[,1]], di[ends_mat[,2]], na.rm = TRUE), # max node distance
                   sum = di[ends_mat[,1]] + di[ends_mat[,2]], # sum node distances (unlikely)
                   harmonic = 2 / (1/di[ends_mat[,1]] + 1/di[ends_mat[,2]]) # weight towards min
  )
  
  d_edge[is.infinite(d_edge)] <- NA
  
  E(g)$edge_dist <- d_edge # estimated prob of edge observation
  
  # Store seed ID 
  graph_attr(g, "seed") <- V(g)$NodeID[seed]
  
  return(g)
}

test1 <- spotlight(Ac2Sim[[1]], seed = 19, combiner = "mean") # function actually works
test2 <- spotlight(Ac2Sim[[1]], seed = 19, combiner = "min")
test3 <- spotlight(Ac2Sim[[1]], seed = 19, combiner = "max")
test4 <- spotlight(Ac2Sim[[1]], seed = 19, combiner = "sum")
test5 <- spotlight(Ac2Sim[[1]], seed = 19, combiner = "harmonic")

V(test)$dist_seed # distances are being accurately stored
V(test)$NodeID

E(test1)$edge_dist
E(test2)$edge_dist
E(test3)$edge_dist
E(test4)$edge_dist
E(test5)$edge_dist

diameter(test) # inverse geodesic reflects actual path lengths

plot(test)

##### Version 3 #####

# Make it so it works on lists of igraph objects, instead of just 1

