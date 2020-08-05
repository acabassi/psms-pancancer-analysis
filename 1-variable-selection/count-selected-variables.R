############################## Selected variables ##############################

rm(list=ls())

for(layer in c("CN", "methylation", "mRNA", "RPPA", "miRNA")){
  for(var_sel in c("_alpha01", "_alpha05", "_alpha1")){
    
    data <- read.csv(paste0("selected-variables/", layer, var_sel, ".csv"),
                     row.names=1)
    cat(layer, var_sel, dim(data)[2], "\n")
  }
  data <- read.csv(paste0("data/", layer, ".csv"), row.names = 1)
  cat(layer, "full", dim(data)[2], "\n")
}
