##################### Pancancer study (Hoadley et al. 2014) ####################
############################## mRNAseq data kernel #############################

rm(list=ls())
library(glmnet)

####################### Load mRNAseq preprocessed data #########################

load("data/preprocessed-mRNA-data.RData")

############################## Load sample names ###############################

load("data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

################################ Select samples ################################

# To check which part of the sample names we need
rownames(mRNAseq)[which((grepl("BT-A20W", rownames(mRNAseq))))]

# Check whether there are any duplicated sample names
sum(duplicated(substr(rownames(mRNAseq), 1, 12))) # None

mRNAseq <- mRNAseq[-which(duplicated(substr(rownames(mRNAseq), 1, 12))),]

# Take only part of sample names that we need 
rownames(mRNAseq) <- substr(rownames(mRNAseq), 1, 12)

# Select only those that are present in Hoadley et al.'s COCA clustering
which_ones <- rownames(mRNAseq)[rownames(mRNAseq)%in%samples]

mRNAseq <- mRNAseq[which_ones, ]
mRNAseq <- mRNAseq[,which(colSums(is.na(mRNAseq))==0)]

dim(mRNAseq)

save(mRNAseq, file = "data/mRNA.csv")
#################### Multinomial regression with elastic-net ###################

response <- annotations_COCA[rownames(mRNAseq),]$Tissue

cv_glmnet <- cv.glmnet(mRNAseq, response, family = "multinomial")
