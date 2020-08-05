##################### Pancancer study (Hoadley et al. 2014) ####################
############################### PSM of RPPA data ###############################

rm(list=ls())

library(impute)
library(klic)
library(pheatmap)
library(PReMiuM)

# chain <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

###################### Load copy number preprocessed data ######################

load("../data/preprocessed-RPPA-data.RData")

############################## Load sample names ###############################

load("../data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

################################ Select samples ################################

# To check which part of the sample names we need
rownames(RPPA)[which((grepl("BT-A20W", rownames(RPPA))))]

# Select only those that are present in Hoadley et al.'s COCA clustering
which_ones <- rownames(RPPA)[rownames(RPPA)%in%samples]

RPPA <- RPPA[which_ones, ]
RPPA_dist <- RPPA_dist[which_ones, which_ones]

save(which_ones, file = "../data/names-RPPA.RData")

RPPA <- as.matrix(RPPA[,6:dim(RPPA)[2]])
RPPA <- impute.knn(RPPA)$data

samples_RPPA <- rownames(RPPA)
save(samples_RPPA, file = "../data/samples_RPPA.RData")

load("../data/samples_intersection.RData")
RPPA <- RPPA[samples_intersection,]

write.csv(RPPA, file = "../data/RPPA.csv", row.names = TRUE)

################################## Run PReMiuM #################################
# 
# n_covariates <- dim(data)[2]
# outcome <- annotations_COCA[rownames(RPPA),]$Tissue
# data <- cbind(outcome, as.data.frame(RPPA))
# 
# prof_regr <-profRegr(yModel="Categorical",
#                      xModel="Normal",
#                      nSweeps=10000,
#                      nClusInit=30,
#                      nBurn=10000,
#                      data=data,
#                      output=paste0("../premium/RPPA_",chain,
#                                    "exclude_y_var_sel"),
#                      covNames = colnames(RPPA),
#                      excludeY = TRUE,
#                      varSelectType = "Continuous",
#                      seed=chain)
# 
# save(prof_regr,
#      file = paste0("../premium/RPPA_chain",
#                    chain,"output_prof_regr_exclude_y_var_sel.RData"))
# 
# dissimObj = PReMiuM::calcDissimilarityMatrix(prof_regr)
# dissMat = PReMiuM::vec2mat(dissimObj$disSimMat, nrow = length(outcome))
# psm <- 1-dissMat
# coph_corr <- copheneticCorrelation(psm)
# 
# save(psm, coph_corr,
#      file = paste0(
#        "../premium/RPPA_chain", chain,"_psm_exclude_y_var_sel.RData"))

# myBlues <-
#     colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)
# png("figures/cc-RPPA.png", height = 800, width = 800)
# pheatmap(cc, annotations = annotations_COCA, color = myBlues, fontsize = 20)
# dev.off()
