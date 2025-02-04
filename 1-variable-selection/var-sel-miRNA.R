##################### Pancancer study (Hoadley et al. 2014) ####################
####################### Variable selection for mRNA data #######################

rm(list=ls())
library(glmnet)

########################## Load mRNA preprocessed data #########################

miRNA <- as.matrix(read.csv("../data/miRNA.csv", header = TRUE, row.names = 1))

############################## Load sample names ###############################

load("../data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

#################### Multinomial regression with elastic-net ###################

response <- factor(annotations_COCA[rownames(miRNA),]$Tissue)
names(response) <- rownames(miRNA)

cv_glmnet_alpha1 <- cv.glmnet(x = miRNA,
                              y = response,
                              family = "multinomial",
                              alpha = 1,
                              standardize = TRUE,
                              type.multinomial = "grouped")

coefficients <- coef(cv_glmnet_alpha1, s = "lambda.min")$BRCA
selected_alpha1 <- rownames(coefficients)[which(coefficients!=0)]

cv_glmnet_alpha05 <- cv.glmnet(x = miRNA,
                               y = response,
                               family = "multinomial",
                               alpha = 0.5,
                               standardize = TRUE,
                               type.multinomial = "grouped")

coefficients <- coef(cv_glmnet_alpha05, s = "lambda.min")$BRCA
selected_alpha05 <-  rownames(coefficients)[which(coefficients!=0)]

cv_glmnet_alpha01 <- cv.glmnet(x = miRNA,
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
     file = "../selected-variables/miRNA.RData")

miRNA_alpha1 <- miRNA[, selected_alpha1[2:51]]
write.csv(miRNA_alpha1, file = "../selected-variables/miRNA_alpha1.RData")
miRNA_alpha05 <- miRNA[, selected_alpha05[2:52]]
write.csv(miRNA_alpha05, file = "../selected-variables/miRNA_alpha05.RData")
miRNA_alpha01 <- miRNA[, selected_alpha01[2:52]]
write.csv(miRNA_alpha01, file = "../selected-variables/miRNA_alpha01.RData")
