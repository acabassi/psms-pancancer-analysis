##################### Pancancer study (Hoadley et al. 2014) ####################
####################### Variable selection for mRNA data #######################

rm(list=ls())
library(glmnet)

########################## Load mRNA preprocessed data #########################

CN <- as.matrix(read.csv("../data/CN.csv", header = TRUE, row.names = 1))

############################## Load sample names ###############################

load("../data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

#################### Multinomial regression with elastic-net ###################

response <- factor(annotations_COCA[rownames(CN),]$Tissue)
names(response) <- rownames(CN)

cv_glmnet_alpha1 <- cv.glmnet(x = CN,
                              y = response,
                              family = "multinomial",
                              alpha = 1,
                              standardize = TRUE,
                              type.multinomial = "grouped")

coefficients <- coef(cv_glmnet_alpha1, s = "lambda.min")$BRCA
selected_alpha1 <- rownames(coefficients)[which(coefficients!=0)]

cv_glmnet_alpha05 <- cv.glmnet(x = CN,
                               y = response,
                               family = "multinomial",
                               alpha = 0.5,
                               standardize = TRUE,
                               type.multinomial = "grouped")

coefficients <- coef(cv_glmnet_alpha05, s = "lambda.min")$BRCA
selected_alpha05 <-  rownames(coefficients)[which(coefficients!=0)]

cv_glmnet_alpha01 <- cv.glmnet(x = CN,
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
     file = "../selected-variables/CN.RData")

CN_alpha1 <- CN[, selected_alpha1[2:85]]
write.csv(CN_alpha1, file = "../selected-variables/CN_alpha1.RData")
CN_alpha05 <- CN[, selected_alpha05[2:85]]
write.csv(CN_alpha05, file = "../selected-variables/CN_alpha05.RData")
CN_alpha01 <- CN[, selected_alpha01[2:85]]
write.csv(CN_alpha01, file = "../selected-variables/CN_alpha01.RData")
