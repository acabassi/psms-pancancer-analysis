### Convergence assessment

rm(list=ls())

# devtools::install_github("acabassi/mdiHelpR")
library(mdiHelpR)
library(stableGR)
library(ggplot2)

n_chains <- 5
noCuda <- "_noCuda"
n_iterations <- 50000

# for(layer in c("CN", "mRNA", "methylation", "miRNA", "RPPA")){ 
layer <- "CN"
  # for(var_sel in c("_alpha1", "_alpha05")){ # , "_alpha01", "_alpha05", "_alpha1"
  var_sel <- "_alpha1"
    
      if(layer %in% c("CN", "mRNA")){
        max_clusters = 250
      }else if(layer %in% c("RPPA", "methylation")){
        max_clusters = 50
      }else if(layer == "miRNA"){
        max_clusters = 100
      }
    
    if(noCuda == "_noCuda" & layer == "mRNA"){
      noCuda <- "_temp"
    }
      
    load(paste0("../mdi/psm_", layer, "_chain", 1, var_sel, noCuda, ".RData"))
    psm_cn1 <- psm_cn; rm(psm_cn)
    load(paste0("../mdi/psm_", layer, "_chain", 2, var_sel, noCuda, ".RData"))
    psm_cn2 <- psm_cn; rm(psm_cn)
    load(paste0("../mdi/psm_", layer, "_chain", 3, var_sel, noCuda, ".RData"))
    psm_cn3 <- psm_cn; rm(psm_cn)
    load(paste0("../mdi/psm_", layer, "_chain", 4, var_sel, noCuda, ".RData"))
    psm_cn4 <- psm_cn; rm(psm_cn)
    load(paste0("../mdi/psm_", layer, "_chain", 5, var_sel, noCuda, ".RData"))
    psm_cn5 <- psm_cn; rm(psm_cn)
      
      # This does not work because there are too many observations
      # compareSimilarityMatrices(psm_cn1,
      #                           psm_cn2,
      #                           psm_cn3,
      #                           psm_cn4,
      #                           psm_cn5)
      # ggsave(paste0("../figures/PSMs_", layer, var_sel,".pdf"))
      
      # psms <- list(psm_cn1, psm_cn2, psm_cn3, psm_cn4, psm_cn5)
      # average_psm <- averageMatrix(psms)
      
      # R hat statistic, Knudson-Vats version
      mcmc_samples_1 <- as.matrix(
        read.csv(paste0("../mdi/", layer, "_", max_clusters, "max_chain", 1,
                        var_sel, "_", n_iterations/1000,"k.csv")
                 )[(n_iterations/10+1):(n_iterations/5),]$MassParameter_1)
      mcmc_samples_2 <- as.matrix(
        read.csv(paste0("../mdi/", layer, "_", max_clusters, "max_chain", 2,
                        var_sel, "_", n_iterations/1000,"k.csv")
                 )[(n_iterations/10+1):(n_iterations/5),]$MassParameter_1)
      mcmc_samples_3 <- as.matrix(
        read.csv(paste0("../mdi/", layer, "_", max_clusters, "max_chain", 3,
                        var_sel, "_", n_iterations/1000,"k.csv")
                 )[(n_iterations/10+1):(n_iterations/5),]$MassParameter_1)
      mcmc_samples_4 <- as.matrix(
        read.csv(paste0("../mdi/", layer, "_", max_clusters, "max_chain", 4,
                        var_sel, "_", n_iterations/1000,"k.csv")
                 )[(n_iterations/10+1):(n_iterations/5),]$MassParameter_1)
      mcmc_samples_5 <- as.matrix(
        read.csv(paste0("../mdi/", layer, "_", max_clusters, "max_chain", 5,
                        var_sel, "_", n_iterations/1000,"k.csv")
                 )[(n_iterations/10+1):(n_iterations/5),]$MassParameter_1)
      stableGR <- stable.GR(list(mcmc_samples_1, mcmc_samples_2, mcmc_samples_3,
                     mcmc_samples_4, mcmc_samples_5))
      reached_target <- (as.numeric(stableGR$psrf) <
                           target.psrf(1, n_chains, epsilon = .1, alpha = .1)$psrf)
      cat("Dataset:", layer, var_sel,"\n",
          "Target reached:", reached_target, "\n\n")
  
#   }
# }

      