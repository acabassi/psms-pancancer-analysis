##################### Pancancer study (Hoadley et al. 2014) ####################
####################### Variable selection for mRNA data #######################

rm(list=ls())
library(glmnet)

########################## Load mRNA preprocessed data #########################

mRNA <- as.matrix(read.csv("../data/mRNA.csv", header = TRUE, row.names = 1))

############################## Load sample names ###############################

load("../data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

#################### Multinomial regression with elastic-net ###################

response <- factor(annotations_COCA[rownames(mRNA),]$Tissue)
names(response) <- rownames(mRNA)

cv_glmnet_alpha1 <- cv.glmnet(x = mRNA,
                              y = response,
                              family = "multinomial",
                              alpha = 1,
                              standardize = TRUE,
                              type.multinomial = "grouped")

coefficients <- coef(cv_glmnet_alpha1, s = "lambda.min")$BRCA
selected_alpha1 <- rownames(coefficients)[which(coefficients!=0)]

cv_glmnet_alpha05 <- cv.glmnet(x = mRNA,
                               y = response,
                               family = "multinomial",
                               alpha = 0.5,
                               standardize = TRUE,
                               type.multinomial = "grouped")

coefficients <- coef(cv_glmnet_alpha05, s = "lambda.min")$BRCA
selected_alpha05 <-  rownames(coefficients)[which(coefficients!=0)]

cv_glmnet_alpha01 <- cv.glmnet(x = mRNA,
                               y = response,
                               family = "multinomial",
                               alpha = 0.1,
                               standardize = TRUE,
                               type.multinomial = "grouped")

coefficients <- coef(cv_glmnet_alpha01, s = "lambda.min")$BRCA
selected_alpha01 <- rownames(coefficients)[which(coefficients!=0)]

save(cv_glmnet_alpha1, selected_alpha1, 
     cv_glmnet_alpha05, selected_alpha05,
     cv_glmnet_alpha01, selected_alpha01,
     file = "../selected-variables/mRNA.RData")

mRNA_alpha1 <- mRNA[, selected_alpha1[2:259]]
write.csv(mRNA_alpha1, file = "../selected-variables/mRNA_alpha1.RData")
mRNA_alpha05 <- mRNA[, selected_alpha05[2:569]]
write.csv(mRNA_alpha05, file = "../selected-variables/mRNA_alpha05.RData")
mRNA_alpha01 <- mRNA[, selected_alpha01[2:1894]]
write.csv(mRNA_alpha01, file = "../selected-variables/mRNA_alpha01.RData")
