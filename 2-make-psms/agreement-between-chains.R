##################### Pancancer study (Hoadley et al. 2014) ####################
###################### Comparison of different MCMC chains #####################

rm(list=ls())

dataset <- "miRNA"
max_clusters <- 100
n_samples <- 50000
noCuda <- "_noCuda" # Can be either "_noCuda" ("_temp" for mRNA) or ""
var_sel <- "_alpha1" 

clusters_list <- list()
for(chain in 1:5){
    load(paste0("../mdi/psm_", dataset, "_chain", chain, var_sel, noCuda,
              ".RData"))
  clusters_list[[chain]] <- clusters
}

ari <- matrix(NA, 5, 5)
diag(ari) <- 1
for(chain_i in 1:4){
  for(chain_j in (chain_i+1):5){
    ari[chain_i,chain_j] <-
      mclust::adjustedRandIndex(clusters_list[[chain_i]],
                                clusters_list[[chain_j]])
  }
}
ari
