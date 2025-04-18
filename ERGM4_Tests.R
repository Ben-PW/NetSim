####################################################################################################

#This file is the test case for all the functions in the ERGM4 version of the code. I'll keep
#them away from the main script to keep things clear. 

####################################################################################################

################################# Randomisation test ###########################################
# Alter this seed to check simulated networks are actually changing
set.seed(2345)

########################################### Begin ##############################################
#### Load environment ####
library(intergraph)
library(dplyr)


#### Importing data ####

source('Data_process.R')
source('Compute_NetStats.R')

####################################### Simulate networks #######################################

FloSim <- simulate(flo1, nsim = 1000)
Ac1Sim <- simulate(acct1, nsim = 1000)
Ac2Sim <- simulate(acct2, nsim = 1000)

rm(flo1, acct1, acct2)
################################ Ensure all networks are undirected ###########################

##### Function to ensure networks are undirected #####

undirect <- function(network_list) {
  lapply(network_list, function(net) {
    if (is.directed(net)) {
      # Convert to adjacency matrix, then back to an undirected network
      net <- as.network(as.matrix.network.adjacency(net), directed = FALSE)
    }
    
    return(net)
  })
}

##### Undirect networks #####

Ac1Sim <- undirect(Ac1Sim)
Ac2Sim <- undirect(Ac2Sim)

###################################### Assign ID to nodes ############################
##### Function to assign node ID #####
assign_node_ids <- function(network_list) {
  for (i in seq_along(network_list)) {
    network_list[[i]] %v% "NodeID" <- seq_len(network.size(network_list[[i]]))
  }
  return(network_list)
}

##### Assign IDs to GT nodes #####
FloSim <- assign_node_ids(FloSim)
Ac1Sim <- assign_node_ids(Ac1Sim)
Ac2Sim <- assign_node_ids(Ac2Sim)

############################### Compute Ground Truth metrics ####################################

##### Function to compute GT metrics for each simulation #####
compute_metrics <- function(network_list, name) {
  data.frame(
    id = paste0(name, seq_along(network_list)),
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

##### Compute GT metrics and combine #####

GTFlo <- compute_metrics(FloSim, "FloSim")
GTAc1 <- compute_metrics(Ac1Sim, "Ac1Sim")
GTAc2 <- compute_metrics(Ac2Sim, "Ac2Sim")

# Combine all into a single dataframe
GTAll <- bind_rows(GTFlo, GTAc1, GTAc2)

# Confirm statistic calculation and storage
cat("Ground truth network statistics calculated and stored in 'GTdf'.\n")

################################### Calculate node level metrics ################################

##### Function to compute and assign node level metrics #####
compute_centrality <- function(network_list) {
  convert_graph <- list()
  
  for(i in seq_along(network_list)) {
    convert_graph[[i]] <- intergraph::asIgraph(network_list[[i]])
    network_list[[i]] %v% "Degree" <- sna::degree(network_list[[i]])
    network_list[[i]] %v% "Betweenness" <- estimate_betweenness(convert_graph[[i]], cutoff = -1)
    network_list[[i]] %v% "Closeness" <- estimate_closeness(convert_graph[[i]], cutoff = -1)
    network_list[[i]] %v% "Eigenvector" <- eigen_centrality(convert_graph[[i]])$vector
    network_list[[i]] %v% "PageRank" <- page_rank(convert_graph[[i]])$vector
  }
  
  return(network_list)
}

##### Compute and assign node level metrics #####

library(igraph)

FloSim <- compute_centrality(FloSim)
Ac1Sim <- compute_centrality(Ac1Sim)
Ac2Sim <- compute_centrality(Ac2Sim)

detach(package:igraph)

######################################### Alter networks ########################################

source("Error_functions_SNA.R") 

###################################### Tie Addition tests #######################################

################## Check correct proportion of ties is actually being added ###################

##### TieAddRand4 #####
# this is the only function that is guaranteed to work every time. Poor performance
# due to manually computing all possible ties. 

# 10% random tie addition
EFlo.01 <- TieAddRand4(FloSim, 0.1)
EAc1.01 <- TieAddRand4(Ac1Sim, 0.1)
EAc2.01 <- TieAddRand4(Ac2Sim, 0.1)

EFlo.01.df <- compute_metrics(EFlo.01, "EFlo01")
EAc1.01.df <- compute_metrics(EAc1.01, "EAc101")
EAc2.01.df <- compute_metrics(EAc2.01, "EAc201")

# 20% random tie addition
EFlo.02 <- TieAddRand4(FloSim, 0.2)
EAc1.02 <- TieAddRand4(Ac1Sim, 0.2) # took 12 seconds
EAc2.02 <- TieAddRand4(Ac2Sim, 0.2) # took 8 seconds

EFlo.02.df <- compute_metrics(EFlo.02,"EFlo01")
EAc1.02df <- compute_metrics(EAc1.02,"EAc101")
EAc2.02df <- compute_metrics(EAc2.02,"EAc201")

# The tie count in these dfs should be 10% and 20% higher than GT

# Updated the TieAddRand4 function to use round() instead of ceiling()
# 10% group 
Flotietest <- (EFlo.01.df$edges - GTFlo$edges)/GTFlo$edges
hist(Flotietest) # all between 7% and 13% increase
Ac1tietest <- (EAc1.01.df$edges - GTAc1$edges)/GTAc1$edges
hist(Ac1tietest) # all between 9.7% and 10.3%
Ac2tietest <- (EAc2.01.df$edges - GTAc2$edges)/GTAc2$edges
hist(Ac2tietest) # all between 9.6% and 10.4%

# 20% group
Flotietest2 <- (EFlo.02.df$edges - GTFlo$edges)/GTFlo$edges
hist(Flotietest2) # all between 16 and 24% increase
Ac1tietest2 <- (EAc1.02df$edges - GTAc1$edges)/GTAc1$edges
hist(Ac1tietest2) # all between 19.7% and 20.2% increase
Ac2tietest2 <- (EAc2.02df$edges - GTAc2$edges)/GTAc2$edges
hist(Ac2tietest2) # all between 19.7% and 20.3% increase

##### TieAddRand6 ######
# This one should be used for larger networks due to use of lookup table. Not 100% guaranteed
# to work but should be fine for anything except huge and incredibly dense networks

# 10% random tie addition
EFlo.01 <- TieAddRand6(FloSim, 0.1)
EAc1.01 <- TieAddRand6(Ac1Sim, 0.1)
EAc2.01 <- TieAddRand6(Ac2Sim, 0.1)

EFlo.01.df <- compute_metrics(EFlo.01, "EFlo01")
EAc1.01.df <- compute_metrics(EAc1.01, "EAc101")
EAc2.01.df <- compute_metrics(EAc2.01, "EAc201")

# 20% random tie addition
EFlo.02 <- TieAddRand6(FloSim, 0.2)
EAc1.02 <- TieAddRand6(Ac1Sim, 0.2) 
EAc2.02 <- TieAddRand6(Ac2Sim, 0.2) 

EFlo.02.df <- compute_metrics(EFlo.02,"EFlo01")
EAc1.02df <- compute_metrics(EAc1.02,"EAc101")
EAc2.02df <- compute_metrics(EAc2.02,"EAc201")

# The tie count in these dfs should be 10% and 20% higher than GT

# Updated the TieAddRand6 function to use round() instead of ceiling()
# 10% group 
Flotietest <- (EFlo.01.df$edges - GTFlo$edges)/GTFlo$edges
hist(Flotietest) # all between 7% and 13% increase
Ac1tietest <- (EAc1.01.df$edges - GTAc1$edges)/GTAc1$edges
hist(Ac1tietest) # all between 9.7% and 10.3%
Ac2tietest <- (EAc2.01.df$edges - GTAc2$edges)/GTAc2$edges
hist(Ac2tietest) # all between 9.6% and 10.4%

# 20% group
Flotietest2 <- (EFlo.02.df$edges - GTFlo$edges)/GTFlo$edges
hist(Flotietest2) # all between 16 and 24% increase
Ac1tietest2 <- (EAc1.02df$edges - GTAc1$edges)/GTAc1$edges
hist(Ac1tietest2) # all between 19.7% and 20.3% increase
Ac2tietest2 <- (EAc2.02df$edges - GTAc2$edges)/GTAc2$edges
hist(Ac2tietest2) # all between 19.7% and 20.3% increase

##### TieAddIterative
# This one will need to be tested to higher proportions to check for adding repeat ties, as I'm
# not sure if that is a risk with this

# 10% random tie addition
EFlo.01 <- TieAddIterative(FloSim, 0.1)
EAc1.01 <- TieAddIterative(Ac1Sim, 0.1)
EAc2.01 <- TieAddIterative(Ac2Sim, 0.1)

EFlo.01.df <- compute_metrics(EFlo.01, "EFlo01")
EAc1.01.df <- compute_metrics(EAc1.01, "EAc101")
EAc2.01.df <- compute_metrics(EAc2.01, "EAc201")

# 20% random tie addition
EFlo.02 <- TieAddIterative(FloSim, 0.2)
EAc1.02 <- TieAddIterative(Ac1Sim, 0.2) 
EAc2.02 <- TieAddIterative(Ac2Sim, 0.2) 

EFlo.02.df <- compute_metrics(EFlo.02,"EFlo01")
EAc1.02df <- compute_metrics(EAc1.02,"EAc101")
EAc2.02df <- compute_metrics(EAc2.02,"EAc201")

# 30% random tie addition
EFlo.03 <- TieAddIterative(FloSim, 0.3)
EAc1.03 <- TieAddIterative(Ac1Sim, 0.3) 
EAc2.03 <- TieAddIterative(Ac2Sim, 0.3) 

EFlo.03.df <- compute_metrics(EFlo.03,"EFlo01")
EAc1.03df <- compute_metrics(EAc1.03,"EAc101")
EAc2.03df <- compute_metrics(EAc2.03,"EAc201")

# 40% random tie addition
EFlo.04 <- TieAddIterative(FloSim, 0.4)
EAc1.04 <- TieAddIterative(Ac1Sim, 0.4) 
EAc2.04 <- TieAddIterative(Ac2Sim, 0.4) 

EFlo.04.df <- compute_metrics(EFlo.04,"EFlo01")
EAc1.04df <- compute_metrics(EAc1.04,"EAc101")
EAc2.04df <- compute_metrics(EAc2.04,"EAc201")

# The tie count in these dfs should be 10% and 20% higher than GT

# Updated the TieAddIterative function to use round() instead of ceiling()
# 10% group 
Flotietest <- (EFlo.01.df$edges - GTFlo$edges)/GTFlo$edges
hist(Flotietest) # all between 7% and 13.5% increase
Ac1tietest <- (EAc1.01.df$edges - GTAc1$edges)/GTAc1$edges
hist(Ac1tietest) # all between 9.7% and 10.3%
Ac2tietest <- (EAc2.01.df$edges - GTAc2$edges)/GTAc2$edges
hist(Ac2tietest) # all between 9.6% and 10.4%

# 20% group
Flotietest2 <- (EFlo.02.df$edges - GTFlo$edges)/GTFlo$edges
hist(Flotietest2) # all between 16 and 24% increase
Ac1tietest2 <- (EAc1.02df$edges - GTAc1$edges)/GTAc1$edges
hist(Ac1tietest2) # all between 19.7% and 20.3% increase
Ac2tietest2 <- (EAc2.02df$edges - GTAc2$edges)/GTAc2$edges
hist(Ac2tietest2) # all between 19.7% and 20.3% increase

# 30% group
Flotietest3 <- (EFlo.03.df$edges - GTFlo$edges)/GTFlo$edges
hist(Flotietest3) # all between 26.5 and 33.5% increase
Ac1tietest3 <- (EAc1.03df$edges - GTAc1$edges)/GTAc1$edges
hist(Ac1tietest3) # all between 29.7% and 30.3% increase
Ac2tietest3 <- (EAc2.03df$edges - GTAc2$edges)/GTAc2$edges
hist(Ac2tietest3) # all between 29.65% and 30.35% increase

# 40% group
Flotietest4 <- (EFlo.04.df$edges - GTFlo$edges)/GTFlo$edges
hist(Flotietest4) # all between 36 and 44% increase
Ac1tietest4 <- (EAc1.04df$edges - GTAc1$edges)/GTAc1$edges
hist(Ac1tietest4) # all between 39.75% and 40.25% increase
Ac2tietest4 <- (EAc2.04df$edges - GTAc2$edges)/GTAc2$edges
hist(Ac2tietest4) # all between 39.7% and 40.3% increase

# Seems like there's a lower bound for function accuracy re size and density

###################################### Tie Removal Tests ######################################

# 10% random tie addition
EFlo.01 <- TieMissRand2(FloSim, 0.1)
EAc1.01 <- TieMissRand2(Ac1Sim, 0.1)
EAc2.01 <- TieMissRand2(Ac2Sim, 0.1)

EFlo.01.df <- compute_metrics(EFlo.01, "EFlo01")
EAc1.01.df <- compute_metrics(EAc1.01, "EAc101")
EAc2.01.df <- compute_metrics(EAc2.01, "EAc201")

# 20% random tie addition
EFlo.02 <- TieMissRand2(FloSim, 0.2)
EAc1.02 <- TieMissRand2(Ac1Sim, 0.2) 
EAc2.02 <- TieMissRand2(Ac2Sim, 0.2) 

EFlo.02.df <- compute_metrics(EFlo.02,"EFlo01")
EAc1.02df <- compute_metrics(EAc1.02,"EAc101")
EAc2.02df <- compute_metrics(EAc2.02,"EAc201")

# 30% random tie addition
EFlo.03 <- TieMissRand2(FloSim, 0.3)
EAc1.03 <- TieMissRand2(Ac1Sim, 0.3) 
EAc2.03 <- TieMissRand2(Ac2Sim, 0.3) 

EFlo.03.df <- compute_metrics(EFlo.03,"EFlo01")
EAc1.03df <- compute_metrics(EAc1.03,"EAc101")
EAc2.03df <- compute_metrics(EAc2.03,"EAc201")

# 40% random tie addition
EFlo.04 <- TieMissRand2(FloSim, 0.4)
EAc1.04 <- TieMissRand2(Ac1Sim, 0.4) 
EAc2.04 <- TieMissRand2(Ac2Sim, 0.4) 

EFlo.04.df <- compute_metrics(EFlo.04,"EFlo01")
EAc1.04df <- compute_metrics(EAc1.04,"EAc101")
EAc2.04df <- compute_metrics(EAc2.04,"EAc201")

# The tie count in these dfs should be 10% and 20% higher than GT

# Updated the TieAddIterative function to use round() instead of ceiling()
par(mfcol = c(3, 4))
# 10% group 
Flotietest <- (EFlo.01.df$edges - GTFlo$edges)/GTFlo$edges
hist(Flotietest) # all between 7% and 13.5% increase
Ac1tietest <- (EAc1.01.df$edges - GTAc1$edges)/GTAc1$edges
hist(Ac1tietest) # all between 9.7% and 10.3%
Ac2tietest <- (EAc2.01.df$edges - GTAc2$edges)/GTAc2$edges
hist(Ac2tietest) # all between 9.6% and 10.4%

# 20% group
Flotietest2 <- (EFlo.02.df$edges - GTFlo$edges)/GTFlo$edges
hist(Flotietest2) # all between 16 and 24% increase
Ac1tietest2 <- (EAc1.02df$edges - GTAc1$edges)/GTAc1$edges
hist(Ac1tietest2) # all between 19.7% and 20.3% increase
Ac2tietest2 <- (EAc2.02df$edges - GTAc2$edges)/GTAc2$edges
hist(Ac2tietest2) # all between 19.7% and 20.3% increase

# 30% group
Flotietest3 <- (EFlo.03.df$edges - GTFlo$edges)/GTFlo$edges
hist(Flotietest3) # all between 26.5 and 33.5% increase
Ac1tietest3 <- (EAc1.03df$edges - GTAc1$edges)/GTAc1$edges
hist(Ac1tietest3) # all between 29.7% and 30.3% increase
Ac2tietest3 <- (EAc2.03df$edges - GTAc2$edges)/GTAc2$edges
hist(Ac2tietest3) # all between 29.65% and 30.35% increase

# 40% group
Flotietest4 <- (EFlo.04.df$edges - GTFlo$edges)/GTFlo$edges
hist(Flotietest4) # all between 36 and 44% increase
Ac1tietest4 <- (EAc1.04df$edges - GTAc1$edges)/GTAc1$edges
hist(Ac1tietest4) # all between 39.75% and 40.25% increase
Ac2tietest4 <- (EAc2.04df$edges - GTAc2$edges)/GTAc2$edges
hist(Ac2tietest4) # all between 39.7% and 40.3% increase

# Testing whether the odd patterns are due to commonalities in simulation
par(mfrow = c(3,1))
hist(GTFlo$edges)
hist(GTAc1$edges)
hist(GTAc2$edges)
# They aren't, honestly I'm happy to continue with them as they are

####################################### Node Removal Tests #################################

# 10% random node deletion
EFlo.01 <- NodeMissRand3(FloSim, 0.1)
EAc1.01 <- NodeMissRand3(Ac1Sim, 0.1)
EAc2.01 <- NodeMissRand3(Ac2Sim, 0.1)

EFlo.01.df <- compute_metrics(EFlo.01, "EFlo01")
EAc1.01.df <- compute_metrics(EAc1.01, "EAc101")
EAc2.01.df <- compute_metrics(EAc2.01, "EAc201")

# 20% random node deletion
EFlo.02 <- NodeMissRand3(FloSim, 0.2)
EAc1.02 <- NodeMissRand3(Ac1Sim, 0.2) 
EAc2.02 <- NodeMissRand3(Ac2Sim, 0.2) 

EFlo.02.df <- compute_metrics(EFlo.02,"EFlo01")
EAc1.02df <- compute_metrics(EAc1.02,"EAc101")
EAc2.02df <- compute_metrics(EAc2.02,"EAc201")

# 30% random node deletion
EFlo.03 <- NodeMissRand3(FloSim, 0.3)
EAc1.03 <- NodeMissRand3(Ac1Sim, 0.3) 
EAc2.03 <- NodeMissRand3(Ac2Sim, 0.3) 

EFlo.03.df <- compute_metrics(EFlo.03,"EFlo01")
EAc1.03df <- compute_metrics(EAc1.03,"EAc101")
EAc2.03df <- compute_metrics(EAc2.03,"EAc201")

# 40% random node deletion
EFlo.04 <- NodeMissRand3(FloSim, 0.4)
EAc1.04 <- NodeMissRand3(Ac1Sim, 0.4) 
EAc2.04 <- NodeMissRand3(Ac2Sim, 0.4) 

EFlo.04.df <- compute_metrics(EFlo.04,"EFlo01")
EAc1.04df <- compute_metrics(EAc1.04,"EAc101")
EAc2.04df <- compute_metrics(EAc2.04,"EAc201")

#histograms wont work because networks retain size after simulation

summary((EFlo.01.df$size - GTFlo$size)/GTFlo$size) #12.5%
summary((EFlo.02.df$size - GTFlo$size)/GTFlo$size) #18.75%
summary((EFlo.03.df$size - GTFlo$size)/GTFlo$size) #31.25%
summary((EFlo.04.df$size - GTFlo$size)/GTFlo$size) #37.5%

#just going to say this function works at this point

#################################### Random Node Addition ####################################

###### Node addition with fixed tie numbers #####

# 10% random node addition
EFlo.01 <- NodeAddFixedEdges2(FloSim, 0.1, 2)
EAc1.01 <- NodeAddFixedEdges2(Ac1Sim, 0.1, 2)
EAc2.01 <- NodeAddFixedEdges2(Ac2Sim, 0.1, 2)

EFlo.01.df <- compute_metrics(EFlo.01, "EFlo01")
EAc1.01.df <- compute_metrics(EAc1.01, "EAc101")
EAc2.01.df <- compute_metrics(EAc2.01, "EAc201")

# 20% random node addition
EFlo.02 <- NodeAddFixedEdges2(FloSim, 0.2, 2)
EAc1.02 <- NodeAddFixedEdges2(Ac1Sim, 0.2, 2) 
EAc2.02 <- NodeAddFixedEdges2(Ac2Sim, 0.2, 2) 

EFlo.02.df <- compute_metrics(EFlo.02,"EFlo01")
EAc1.02df <- compute_metrics(EAc1.02,"EAc101")
EAc2.02df <- compute_metrics(EAc2.02,"EAc201")

par(mfcol = c(3, 2))
# 10% group 
Flotietest <- (EFlo.01.df$edges - GTFlo$edges)/GTFlo$edges
hist(Flotietest) # all between 7% and 13.5% increase
Ac1tietest <- (EAc1.01.df$edges - GTAc1$edges)/GTAc1$edges
hist(Ac1tietest) # all between 9.7% and 10.3%
Ac2tietest <- (EAc2.01.df$edges - GTAc2$edges)/GTAc2$edges
hist(Ac2tietest) # all between 9.6% and 10.4%

# 20% group
Flotietest2 <- (EFlo.02.df$edges - GTFlo$edges)/GTFlo$edges
hist(Flotietest2) # all between 16 and 24% increase
Ac1tietest2 <- (EAc1.02df$edges - GTAc1$edges)/GTAc1$edges
hist(Ac1tietest2) # all between 19.7% and 20.3% increase
Ac2tietest2 <- (EAc2.02df$edges - GTAc2$edges)/GTAc2$edges
hist(Ac2tietest2) # all between 19.7% and 20.3% increase

# The edge data for all the ones except Flo look weird
# Ac series all have approximately 10% of the increase Flo series have
# However this could be due to the fixed number of added ties and relatively small added nodes?
# This makes sense, Flo series will receive approx 1 extra node, typically containing between
# 14 and 24 ties, so the addition of 1 node and two ties is significant
# AC series have 26 nodes, but approx 180 ties, meaning approx 2 nodes added with 4 ties, 
# which is obviously minimal.
# Time to break out the NodeAddRandDD function

###### Add nodes with degree randomly based on degree distribution #####
par(mfcol = c(1,1))
plot(FloSim[[1]])
FloTest <- NodeAddRandDD(FloSim[1], 0.2)
plot(FloTest[[1]])
network::is.multiplex(FloSim[[1]])

network.size(FloTest[[1]])
network.size(FloSim[[1]])
degree(FloSim[[1]], gmode = "graph")
degree(FloTest[[1]], gmode = "graph")
#histograms wont work because networks retain size after simulation

summary((EFlo.01.df$size - GTFlo$size)/GTFlo$size) #12.5%
summary((EFlo.02.df$size - GTFlo$size)/GTFlo$size) #18.75%
summary((EFlo.03.df$size - GTFlo$size)/GTFlo$size) #31.25%
summary((EFlo.04.df$size - GTFlo$size)/GTFlo$size) #37.5%

#just going to say this function works at this point

######################################### Store data ###########################################

##### Global network metrics #####
EFlo.01.df <- compute_metrics(EFlo.01, "EFlo")
EAc1.01.df <- compute_metrics(EAc1.01, "EAc1")
EAc2.01.df <- compute_metrics(EAc2.01, "EAc2")
EAc1NodeAdd.df <- compute_metrics(EAc1NodeAdd, "EAc1NodeAdd")

##### Node level metrics #####
library(igraph)

EFlo.01 <- compute_centrality(EFlo.01)
EAc1.01 <- compute_centrality(EAc1.01)
EAc2.01 <- compute_centrality(EAc2.01)
ENodeAdd <- compute_centrality(EAc1NodeAdd)

detach(package:igraph)


################################# Dataframe to compare centrality scores ##########################

# IMPORTANT: I have removed the code squaring the pearson correlation because I'm not 100% on 
# the justification for squaring it

##### Function to compute all relevant centrality metrics and store in dataframe #####

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

##### Compute bias #####

# Compute bias for all network simulations
bias_FloSim <- compute_bias_metrics7(FloSim, EFlo.01, "FloSim")
bias_FloSim1 <- compute_bias_metrics7(FloSim, EFlo.01, "FloSimTest")
bias_Ac1Sim <- compute_bias_metrics7(Ac1Sim, EAc1.01, "Ac1Sim")
bias_Ac2Sim <- compute_bias_metrics7(Ac2Sim, EAc2.01, "Ac2Sim")
bias_NodeAdd <- compute_bias_metrics7(Ac1Sim, ENodeAdd, "NodeTest")


# Combine all results into a single dataframe
metric_bias <- bind_rows(bias_FloSim, bias_Ac1Sim, bias_Ac2Sim, bias_NodeAdd)

# Add in original full network ground truth metrics
bias_GTadd <- left_join(metric_bias, GTAll, by = "id")

