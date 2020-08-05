##################### Pancancer study (Hoadley et al. 2014) ####################
#################### Variable selection for methylation data ###################

rm(list=ls())
library(glmnet)

###################### Load methylation preprocessed data ######################

methylation <-
  as.matrix(read.csv("../data/methylation.csv", header = TRUE, row.names = 1))

############################## Load sample names ###############################

load("../data/samples.RData")
annotations_COCA <- anno_col; rm(anno_col)

################### Multinomial regression with elastic-net ###################

response <- factor(annotations_COCA[rownames(methylation),]$Tissue)
names(response) <- rownames(methylation)

cv_glmnet_alpha1 <- cv.glmnet(x = methylation,
                              y = response,
                              family = "multinomial",
                              alpha = 1,
                              standardize = FALSE,
                              type.multinomial = "grouped")

plot(cv_glmnet_alpha1)
cv_glmnet_alpha1$lambda.min

coefficients <- coef(cv_glmnet_alpha1, s = "lambda.min")$BRCA
selected_alpha1 <- rownames(coefficients)[which(coefficients!=0)]
length(selected_alpha1)

cv_glmnet_alpha05 <- cv.glmnet(x = methylation,
                              y = response,
                              family = "multinomial",
                              alpha = 0.5,
                              standardize = FALSE,
                              type.multinomial = "grouped")

plot(cv_glmnet_alpha05)
cv_glmnet_alpha05$lambda.min

coefficients <- coef(cv_glmnet_alpha05, s = "lambda.min")$BRCA
selected_alpha05 <- rownames(coefficients)[which(coefficients!=0)]

cv_glmnet_alpha01 <- cv.glmnet(x = methylation,
                               y = response,
                               family = "multinomial",
                               alpha = 0.1,
                               standardize = FALSE,
                               type.multinomial = "grouped")

plot(cv_glmnet_alpha01)
cv_glmnet_alpha01$lambda.min

coefficients <- coef(cv_glmnet_alpha01, s = "lambda.min")$BRCA
selected_alpha01 <- rownames(coefficients)[which(coefficients!=0)]

save(cv_glmnet_alpha1, cv_glmnet_alpha05, cv_glmnet_alpha01,
     selected_alpha1, selected_alpha05, selected_alpha01,
     file = "../selected-variables/methylation.RData")

methylation_alpha1 <- methylation[, selected_alpha1[2:323]]
write.csv(methylation_alpha1,
          file = "../selected-variables/methylation_alpha1.RData")
methylation_alpha05 <- methylation[, selected_alpha05[2:624]]
write.csv(methylation_alpha05,
          file = "../selected-variables/methylation_alpha05.RData")
methylation_alpha01 <- methylation[, selected_alpha01[2:1440]]
write.csv(methylation_alpha01,
          file = "../selected-variables/methylation_alpha01.RData")
