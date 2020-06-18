##################### Pancancer study (Hoadley et al. 2014) ####################
############################ Copy number data kernel ###########################

rm(list=ls())

library(impute)
library(klic)
library(pheatmap)
library(PReMiuM)

# chain <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

###################### Load copy number preprocessed data ######################

load("../data/preprocessed-SCNA-data.RData")

############################## Load sample names ###############################

load("../data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

################################ Select samples ################################

# To check which part of the sample names we need
rownames(SCNAscaled)[which((grepl("BT-A20W", rownames(SCNAscaled))))]

# Check whether there are any duplicated sample names
which(duplicated(substr(rownames(SCNAscaled), 1, 12))) # None

# Take only part of sample names that we need 
rownames(SCNAscaled) <- substr(rownames(SCNAscaled), 1, 12)

# Select only those that are present in Hoadley et al.'s COCA clustering
which_ones <- rownames(SCNAscaled)[rownames(SCNAscaled)%in%samples]

SCNAscaled <- SCNAscaled[which_ones, ]
save(which_ones, file = "../data/names-CN.RData")

colnames(SCNAscaled)[42] <- "4p16.3.2"
colnames(SCNAscaled) <- make.names(colnames(SCNAscaled))
# miceOut <- mice(SCNAscaled)
# SCNAscaled <- complete(miceOut)
SCNAscaled <- impute.knn(SCNAscaled)$data

samples_CN <- rownames(SCNAscaled)
save(samples_CN, file = "../data/samples_CN.RData")

load("../data/samples_intersection.RData")

SCNAscaled <- SCNAscaled[samples_intersection,]

write.csv(SCNAscaled, file = "../data/CN.csv", row.names = TRUE)

################################ Run PReMiuM ###################################

# n_covariates <- dim(data)[2]
# outcome <- annotations_COCA[rownames(SCNAscaled),]$Tissue
# data <- cbind(outcome, as.data.frame(SCNAscaled))
# 
# prof_regr <-profRegr(yModel="Categorical",
#                      xModel="Normal",
#                      nSweeps=10000,
#                      nClusInit=30,
#                      nBurn=10000,
#                      data=data,
#                      output=paste0("../premium/CN_",chain,"exclude_y_var_sel"),
#                      covNames = colnames(SCNAscaled),
#                      excludeY = TRUE,
#                      varSelectType = "Continuous",
#                      seed=chain)
# 
# save(prof_regr,
#      file = paste0("../premium/CN_chain",
#                    chain,"output_prof_regr_exclude_y_var_sel.RData"))
# 
# dissimObj = PReMiuM::calcDissimilarityMatrix(prof_regr)
# dissMat = PReMiuM::vec2mat(dissimObj$disSimMat, nrow = length(outcome))
# psm <- 1-dissMat
# coph_corr <- copheneticCorrelation(psm)
# 
# save(psm, coph_corr,
#      file = paste0("../premium/CN_chain", chain,"_psm_exclude_y_var_sel.RData"))

# myBlues <-
#     colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)
# png("figures/cc-CN.png", height = 800, width = 800)
# pheatmap(cc, annotations = annotations_COCA, color = myBlues, fontsize = 20)
# dev.off()
