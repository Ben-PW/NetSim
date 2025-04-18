##### Bias calculation 2 #####

compute_bias_metrics8 <- function(original_sim, perturbed_sim, name) {
  
  for(i in seq_along(original_sim)){
  
  common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
  
  return(common_ids)
  
  }
  
  data.frame(
    id = paste0(name, seq_along(original_sim)),
    
    deg_robust_cor = sapply(seq_along(original_sim), function(i) {
      # Match nodes by NodeID
      #common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor((original_sim[[i]] %v% "Degree")[orig_indices], 
          (perturbed_sim[[i]] %v% "Degree")[pert_indices], 
          use = "complete.obs")
    }),
    
    bet_robust_cor = sapply(seq_along(original_sim), function(i) {
      #common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor((original_sim[[i]] %v% "Betweenness")[orig_indices], 
          (perturbed_sim[[i]] %v% "Betweenness")[pert_indices], 
          use = "complete.obs")
    }),
    
    clo_robust_cor = sapply(seq_along(original_sim), function(i) {
      #common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor((original_sim[[i]] %v% "Closeness")[orig_indices], 
          (perturbed_sim[[i]] %v% "Closeness")[pert_indices], 
          use = "complete.obs")
    }),
    
    eig_robust_cor = sapply(seq_along(original_sim), function(i) {
      #common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor((original_sim[[i]] %v% "Eigenvector")[orig_indices], 
          (perturbed_sim[[i]] %v% "Eigenvector")[pert_indices], 
          use = "complete.obs")
    }),
    
    pg_robust_cor = sapply(seq_along(original_sim), function(i) {
      #common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor((original_sim[[i]] %v% "PageRank")[orig_indices], 
          (perturbed_sim[[i]] %v% "PageRank")[pert_indices], 
          use = "complete.obs")
    }),
    
    deg_robust_srcc = sapply(seq_along(original_sim), function(i) {
      #common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor(rank((original_sim[[i]] %v% "Degree")[orig_indices]), 
          rank((perturbed_sim[[i]] %v% "Degree")[pert_indices]), 
          method = "spearman", use = "complete.obs")
    }),
    
    bet_robust_srcc = sapply(seq_along(original_sim), function(i) {
      #common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor(rank((original_sim[[i]] %v% "Betweenness")[orig_indices]), 
          rank((perturbed_sim[[i]] %v% "Betweenness")[pert_indices]), 
          method = "spearman", use = "complete.obs")
    }),
    
    clo_robust_srcc = sapply(seq_along(original_sim), function(i) {
      #common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor(rank((original_sim[[i]] %v% "Closeness")[orig_indices]), 
          rank((perturbed_sim[[i]] %v% "Closeness")[pert_indices]), 
          method = "spearman", use = "complete.obs")
    }),
    
    eig_robust_srcc = sapply(seq_along(original_sim), function(i) {
      #common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor(rank((original_sim[[i]] %v% "Eigenvector")[orig_indices]), 
          rank((perturbed_sim[[i]] %v% "Eigenvector")[pert_indices]), 
          method = "spearman", use = "complete.obs")
    }),
    
    pg_robust_srcc = sapply(seq_along(original_sim), function(i) {
      #common_ids <- intersect(original_sim[[i]] %v% "NodeID", perturbed_sim[[i]] %v% "NodeID")
      orig_indices <- match(common_ids, original_sim[[i]] %v% "NodeID")
      pert_indices <- match(common_ids, perturbed_sim[[i]] %v% "NodeID")
      
      cor(rank((original_sim[[i]] %v% "PageRank")[orig_indices]), 
          rank((perturbed_sim[[i]] %v% "PageRank")[pert_indices]), 
          method = "spearman", use = "complete.obs")
    })
  )
}
