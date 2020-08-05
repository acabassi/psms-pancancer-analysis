##################### Pancancer study (Hoadley et al. 2014) ####################
####################### Variable selection for mRNA data #######################

rm(list=ls())
library(glmnet)

########################## Load mRNA preprocessed data #########################

RPPA <- as.matrix(read.csv("../data/RPPA.csv", header = TRUE, row.names = 1))

############################## Load sample names ###############################

load("../data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

#################### Multinomial regression with elastic-net ###################

response <- factor(annotations_COCA[rownames(RPPA),]$Tissue)
names(response) <- rownames(RPPA)

cv_glmnet_alpha1 <- cv.glmnet(x = RPPA,
                              y = response,
                              family = "multinomial",
                              alpha = 1,
                              standardize = TRUE,
                              type.multinomial = "grouped")

coefficients <- coef(cv_glmnet_alpha1, s = "lambda.min")$BRCA
selected_alpha1 <- rownames(coefficients)[which(coefficients!=0)]

cv_glmnet_alpha05 <- cv.glmnet(x = RPPA,
                               y = response,
                               family = "multinomial",
                               alpha = 0.5,
                               standardize = TRUE,
                               type.multinomial = "grouped")

coefficients <- coef(cv_glmnet_alpha05, s = "lambda.min")$BRCA
selected_alpha05 <-  rownames(coefficients)[which(coefficients!=0)]

cv_glmnet_alpha01 <- cv.glmnet(x = RPPA,
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
     file = "../selected-variables/RPPA.RData")

RPPA_alpha1 <- RPPA[, selected_alpha1[2:125]]
write.csv(RPPA_alpha1, file = "../selected-variables/RPPA_alpha1.RData")
RPPA_alpha05 <- RPPA[, selected_alpha05[2:132]]
write.csv(RPPA_alpha05, file = "../selected-variables/RPPA_alpha05.RData")
RPPA_alpha01 <- RPPA[, selected_alpha01[2:132]]
write.csv(RPPA_alpha01, file = "../selected-variables/RPPA_alpha01.RData")
