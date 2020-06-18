##################### Pancancer study (Hoadley et al. 2014) ####################
############################ PSM of methylation data ###########################

rm(list=ls())

library(bnstruct)
library(klic)
library(parallel)
# library(pheatmap)
# library(PReMiuM)

# chain <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

###################### Load methylation preprocessed data ######################

load("../data/preprocessed-methylation-data.RData")

############################## Load sample names ###############################

load("../data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

################################ Select samples ################################

# To check which part of the sample names we need
rownames(methylation)[which((grepl("BT-A20W", rownames(methylation))))]

# Check whether there are any duplicated sample names
sum(duplicated(substr(rownames(methylation), 1, 12))) # None

methylation <-
  methylation[-which(duplicated(substr(rownames(methylation), 1, 12))),]
methylation_dist <- as.matrix(methylation_dist)
methylation_dist <- methylation_dist[
  -which(duplicated(substr(rownames(methylation_dist), 1, 12))),
  -which(duplicated(substr(colnames(methylation_dist), 1, 12)))]

# Take only part of sample names that we need 
rownames(methylation) <- substr(rownames(methylation), 1, 12)
rownames(methylation_dist) <- colnames(methylation_dist) <-
  substr(rownames(methylation_dist), 1, 12)

# Select only those that are present in Hoadley et al.'s COCA clustering
which_ones <- rownames(methylation)[rownames(methylation)%in%samples]

methylation <- methylation[which_ones, ]
methylation_dist <- methylation_dist[which_ones, which_ones]

# save(which_ones, file = "../data/names-methylation.RData")
# 
# system.time(imp <-
#               parlmice(data = methylation[,1:200], n.core = 5, n.imp.core = 1,
#                        method = "logreg"))
# 
methylation <- as.data.frame(methylation)

samples_methylation <- rownames(methylation)
save(samples_methylation, file = "../data/samples_methylation.RData")

methylation <- sapply(methylation, as.factor)
methylation <- knn.impute(methylation)
rownames(methylation) <- samples_methylation

load("../data/samples_intersection.RData")

methylation <- methylation[samples_intersection,]

write.csv(methylation, file = "../data/methylation.csv", row.names = TRUE)

################################## Run PReMiuM #################################
# 
# n_covariates <- dim(data)[2]
# outcome <- annotations_COCA[rownames(methylation),]$Tissue
# data <- cbind(outcome, as.data.frame(methylation))
# 
# prof_regr <-profRegr(yModel="Categorical",
#                      xModel="Discrete",
#                      nSweeps=100000,
#                      nClusInit=30,
#                      nBurn=20000,
#                      data=data,
#                      output=paste0("../premium/methylation_",chain,
#                                    "exclude_y_var_sel"),
#                      covNames = colnames(methylation),
#                      excludeY = TRUE,
#                      varSelectType = "Continuous",
#                      seed=chain)
# 
# save(prof_regr,
#      file = paste0("../premium/methylation_chain",
#                    chain,"output_prof_regr_exclude_y_var_sel.RData"))
# 
# dissimObj = PReMiuM::calcDissimilarityMatrix(prof_regr)
# dissMat = PReMiuM::vec2mat(dissimObj$disSimMat, nrow = length(outcome))
# psm <- 1-dissMat
# coph_corr <- copheneticCorrelation(psm)
# 
# save(psm, coph_corr,
#      file = paste0(
#        "../premium/methylation_chain", chain,"_psm_exclude_y_var_sel.RData"))

# myBlues <-
#     colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)
# png("figures/cc-methylation.png", height = 800, width = 800)
# pheatmap(cc, annotations = annotations_COCA, color = myBlues, fontsize = 20)
# dev.off()
