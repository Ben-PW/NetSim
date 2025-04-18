####################################################################################################

#Combine3.R is the updated version of Combine2.R. The aims of Combine3 are to begin compartmentalising
#and updating the code. Data import and ERGM fitting will be performed in separate scripts. 
#Likewise, data visualisation etc will be performed in separate scripts. The main improvement for this
#iteration is to update the code such that it can handle simulated networks from multiple different
#ERGMs. An ID variable referencing the ERGM name and simulation iteration will be necessary

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

##### Simple tie missingness #####
source("Error_functions.R") 

EFlo.01 <- TieMissRand(FloSim, 0.1) # Only 1 [] required if altering 1 network
EAc1.01 <- TieMissRand(Ac1Sim, 0.1)
EAc2.01 <- TieMissRand(Ac2Sim, 0.1)
EFlo.02 <- TieAddRand(FloSim, 0.1)
EFlo.03 <- NodeAddRand_DD(FloSim, 0.1)

##### TEST FOR NODE ADD FUNCTIONALITY #####

detach(package:igraph)
EfloTest <- NodeAddRand_DD(FloSim[1], 0.1)
FloTestdf <- compute_metrics(EFlo.03, "Test")

######################################### Store data ###########################################

##### Global network metrics #####
EFlo.01.df <- compute_metrics(EFlo.01, "EFlo")
EAc1.01.df <- compute_metrics(EAc1.01, "EAc1")
EAc2.01.df <- compute_metrics(EAc2.01, "EAc2")

##### Node level metrics #####
library(igraph)

EFlo.01 <- compute_centrality(EFlo.01)
EAc1.01 <- compute_centrality(EAc1.01)
EAc2.01 <- compute_centrality(EAc2.01)
ENodeAdd <- compute_centrality(EFlo.03)

detach(package:igraph)

##### Cheeky little data check #####
# Degree is easy to check if ties are removed - the error network can be subtracted from the 
# ground truth and there should never be a negative value

check <- list()
for (i in seq_along(FloSim)) {
  check[[i]] <- FloSim[[i]] %v% "Degree" - EFlo.01[[i]] %v% "Degree"
}

if (all(check[[i]] >= 0)) {
  print("We *probably* Gucci")
  rm(check)
} else {
  print("We are decidedly not Gucci")
}
# (they were, in fact, Gucci)
################################# Dataframe to compare centrality scores ##########################

# IMPORTANT: I have removed the code squaring the pearson correlation because I'm not 100% on 
# the justification for squaring it

##### Function to compute all relevant centrality metrics and store in dataframe #####

compute_bias_metrics <- function(original_sim, perturbed_sim, name) {
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

##### Compute bias #####

# Compute bias for all network simulations
bias_FloSim <- compute_bias_metrics(FloSim, EFlo.01, "FloSim")
bias_FloSim1 <- compute_bias_metrics3(FloSim, EFlo.01, "FloSimTest")
bias_Ac1Sim <- compute_bias_metrics(Ac1Sim, EAc1.01, "Ac1Sim")
bias_Ac2Sim <- compute_bias_metrics(Ac2Sim, EAc2.01, "Ac2Sim")
bias_NodeAdd <- compute_bias_metrics3(FloSim, ENodeAdd, "NodeTest")


# Combine all results into a single dataframe
metric_bias <- bind_rows(bias_FloSim, bias_Ac1Sim, bias_Ac2Sim)

# Add in original full network ground truth metrics
bias_GTadd <- left_join(metric_bias, GTAll, by = "id")


########################################### Visualisations ######################################

library(ggplot2)
library(reshape2)

bias_long <- melt(metric_bias, 
                  id.vars = "id", 
                  variable.name = "Centrality", 
                  value.name = "Robustness")

ggplot(data = bias_long,
       aes(x = Centrality,
           y = Robustness)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(bias_GTadd, aes(x = density, y = bet_robust_srcc)) +
  geom_point(alpha = 0.6, color = "blue") +
  theme_minimal() 
 

