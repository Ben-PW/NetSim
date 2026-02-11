#################################################################################################

# Complete re-write of spotlight function

####################################################################################

assignSpotlight <- function(graph_list, spotlight_pct, s) {
  lapply(graph_list, function(g) {
    
    n <- vcount(g)
    k <- round(n * spotlight_pct)
    
    Spotlight <- rep(0, n) # initialise all with no spotlight
    
    deg <- degree(g) # extract node degrees
    w   <- (1 - s) * rep(1, n) + s * deg # weight sampling by degree and s
    
    idx <- sample(seq_len(n), k, prob = w)
    Spotlight[idx] <- 1 # set selected node spotlights to 1
    
    V(g)$Spotlight <- Spotlight
    g
  })
}


