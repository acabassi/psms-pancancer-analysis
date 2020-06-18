#################### Pancancer study (Hoadley et al. 2014) #####################
############################# PSM of miRNAseq data #############################

rm(list=ls())

library(impute)
library(klic)
library(pheatmap)
library(PReMiuM)

# chain <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

####################### Load miRNAseq preprocessed data ########################

load("../data/preprocessed-miRNA-data.RData")
miRNA <- final_miRNA_data; rm(final_miRNA_data)

############################# Load sample names ################################

load("../data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

################################ Select samples ################################

# To check which part of the sample names we need
rownames(miRNA)[which((grepl("BT-A20W", rownames(miRNA))))]

# Check whether there are any duplicated sample names
sum(duplicated(substr(rownames(miRNA), 6, 17))) # None

miRNA <- miRNA[-which(duplicated(substr(rownames(miRNA), 6, 17))),]
miRNA_dist <- as.matrix(miRNA_dist)
miRNA_dist <- miRNA_dist[
  -which(duplicated(substr(rownames(miRNA_dist), 6, 17))), 
  -which(duplicated(substr(rownames(miRNA_dist), 6, 17)))]

# Take only part of sample names that we need 
rownames(miRNA) <- substr(rownames(miRNA), 6, 17)
rownames(miRNA_dist) <- colnames(miRNA_dist) <-
  substr(rownames(miRNA_dist), 6, 17)

# Select only those that are present in Hoadley et al.'s COCA clustering
which_ones <- rownames(miRNA)[rownames(miRNA)%in%samples]

miRNA <- miRNA[which_ones, ]
miRNA_dist <- miRNA_dist[which_ones, which_ones]

save(which_ones, file = "../data/names-miRNA.RData")

miRNA <- impute.knn(miRNA)$data

samples_miRNA <- rownames(miRNA)
save(samples_miRNA, file = "../data/samples_miRNA.RData")

load("../data/samples_intersection.RData")
miRNA <- miRNA[samples_intersection,]

write.csv(miRNA, file = "../data/miRNA.csv", row.names = TRUE)

################################## Run PReMiuM #################################
# 
# n_covariates <- dim(data)[2]
# outcome <- annotations_COCA[rownames(miRNA),]$Tissue
# data <- cbind(outcome, as.data.frame(miRNA))
# 
# prof_regr <-profRegr(yModel="Categorical",
#                      xModel="Normal",
#                      nSweeps=10000,
#                      nClusInit=30,
#                      nBurn=10000,
#                      data=data,
#                      output=paste0("../premium/miRNA_",chain,
#                                    "exclude_y_var_sel"),
#                      covNames = colnames(miRNA),
#                      excludeY = TRUE,
#                      varSelectType = "Continuous",
#                      seed=chain)
# 
# save(prof_regr,
#      file = paste0("../premium/miRNA_chain",
#                    chain,"output_prof_regr_exclude_y_var_sel.RData"))
# 
# dissimObj = PReMiuM::calcDissimilarityMatrix(prof_regr)
# dissMat = PReMiuM::vec2mat(dissimObj$disSimMat, nrow = length(outcome))
# psm <- 1-dissMat
# coph_corr <- copheneticCorrelation(psm)
# 
# save(psm, coph_corr,
#      file = paste0(
#        "../premium/miRNA_chain", chain,"_psm_exclude_y_var_sel.RData"))

# myBlues <-
#     colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)
# png("figures/cc-miRNA.png", height = 800, width = 800)
# pheatmap(cc, annotations = annotations_COCA, color = myBlues, fontsize = 20)
# dev.off()
