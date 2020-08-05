####################### Transcriptional module discovery #######################
########################### Unsupervised integration ###########################

rm(list = ls())

set.seed(1)

library(circlize)
library(ComplexHeatmap)
library(colorspace)
library(mclust)
library(rcartocolor)
library(R.matlab)

# library(devtools)
# install_github("acabassi/coca")
library(coca)
# install_github("acabassi/klic")
library(klic)

var_sel <- "_alpha05"
noCuda <- "_noCuda"
if(noCuda == "_noCuda"){
  noCuda_mRNA = "_temp"
}else{
  noCuda_mRNA = ""
}

################################# Load PSMs  ###################################

load(paste0("../mdi/psms-rdata/psm_CN_average_",
            var_sel, noCuda, ".RData"))
psm_cn <- average_psm; rm(average_psm)
load(paste0("../mdi/psms-rdata/psm_methylation_average_",
            var_sel, noCuda, ".RData"))
psm_methylation <- average_psm; rm(average_psm)
load(paste0("../mdi/psms-rdata/psm_miRNA_average_",
            var_sel, noCuda, ".RData"))
psm_mirna <- average_psm; rm(average_psm)
load(paste0("../mdi/psms-rdata/psm_mRNA_average_",
            var_sel, noCuda_mRNA, ".RData"))
psm_mrna <- average_psm; rm(average_psm)
load(paste0("../mdi/psms-rdata/psm_RPPA_average_",
            var_sel, noCuda, ".RData"))
psm_rppa <- average_psm; rm(average_psm)

############## Check that statistical units are in the same order ##############

sum(1-(rownames(psm_cn) == rownames(psm_methylation))) # Should be zero
sum(1-(rownames(psm_cn) == rownames(psm_mirna))) # Should be zero
sum(1-(rownames(psm_cn) == rownames(psm_mrna))) # Should be zero
sum(1-(rownames(psm_cn) == rownames(psm_rppa))) # Should be zero

################################ Spectrum shift ################################

psm_cn <- spectrumShift(psm_cn, verbose = TRUE, coeff = 1.01)
psm_methylation <- spectrumShift(psm_methylation, verbose = TRUE, coeff = 1.01)
psm_mirna <- spectrumShift(psm_mirna, verbose = TRUE, coeff = 1.01)
psm_mrna <- spectrumShift(psm_mrna, verbose = TRUE, coeff = 1.01)
psm_rppa <- spectrumShift(psm_rppa, verbose = TRUE, coeff = 1.01)

######################### Create array of all matrices #########################

n_datasets <- 5
datasetNames <- c("CN", "Methylation", "miRNA", "mRNA", "RPPA")
N <- dim(psm_cn)[1]
CM <- array(NA, c(N, N, n_datasets))
dimnames(CM) <- list(rownames(psm_cn), rownames(psm_cn), datasetNames)

CM[,, "CN"] <- psm_cn
CM[,, "Methylation"] <- psm_methylation
CM[,, "miRNA"] <- psm_mirna
CM[,, "mRNA"] <- psm_mrna
CM[,, "RPPA"] <- psm_rppa

################################# Integration ##################################
# Use localised multiple kernel k-means to integrate the datasets

# Initialise array of kernel matrices
maxK <- 50
KM <- array(0, c(N, N, maxK - 1))
clLabels <- array(NA, c(maxK - 1, N))
theta <- array(NA, c(maxK - 1, N, n_datasets))

### Each of these must is run separately on the HPC ###
parameters <- list()
parameters$iteration_count <-
    500 # set the maximum number of iterations

# i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#     # Use kernel k-means with K=i to find weights and cluster labels
#     parameters$cluster_count <- i # set the number of clusters K
#     lmkkm <- lmkkmeans(CM, parameters, verbose = TRUE)
# 
#     # Compute weighted matrix
#     for (j in 1:dim(CM)[3]) {
#         KM[, , i - 1] <-
#             KM[, , i - 1] +
#             (lmkkm$Theta[, j] %*% t(lmkkm$Theta[, j])) * CM[, , j]
#     }
# 
#     # Save cluster labels
#     clLabels[i - 1, ] <- lmkkm$clustering
#     theta[i - 1, , ] <- lmkkm$Theta
# 
# KM_i <- KM[,,i-1]
# clLabels_i <- clLabels[i-1,]
# theta_i <- theta[i-1,,]

# save(KM_i, clLabels_i, theta_i,
#      file = paste0("integration-output/unsupervised_output_", i, ".RData"))

for(i in 2:maxK){
    load(paste0("../unsupervised-integration-output/unsupervised_output_",
    i, var_sel, ".RData"))
    KM[,,i-1] <- KM_i
    clLabels[i-1,] <- clLabels_i
    theta[i-1,,] <- theta_i
}
    
############################### Maximise silhouette ############################

maxSil <- maximiseSilhouette(KM, clLabels, maxK = maxK)

save(
    KM,
    clLabels,
    theta,
    maxSil,
    file = paste0(
      "../unsupervised-integration-output/unsupervised_output_all",
      var_sel,".RData")
)

# load("data/unsupervised_output_dpmsysbio.RData")

dimnames(theta) <- list(paste(as.character(2:maxK), "_clusters"),
                        rownames(CM),
                        c("CN", "Methylation", "miRNA", "mRNA", "RPPA")
)

bestK <- maxSil$K
inds <- rownames(CM)
clustersBestK <- clLabels[bestK - 1, ]
names(clustersBestK) <- rownames(CM)

################################## Plot output #################################

# Let's just plot the clusters:
table(clustersBestK)
sortedClusters <- sort(clustersBestK, index.return = T)

myBlues <- colorRamp2(c(0, 1), colors = c("white", "#003C71"))
my_anno_legend_param <- function(name, nrow=1) {
  return(list(
    title = name,
    labels_gp = gpar(fontsize = 18),
    title_gp = gpar(fontsize = 22),
    nrow = nrow
  ))
}

col_bw = c("white", "black")

# This contains palette12, which is used for the tumour colours throughout
load("../data/tumour_colours.RData")

names(clustersBestK) <- gsub("\\.", "-", names(clustersBestK))
load("../data/samples.RData")
adjustedRandIndex(anno_col[names(clustersBestK),]$Tissue, clustersBestK)
adjustedRandIndex(anno_col[names(clustersBestK),]$COCA, clustersBestK)

# I ran this only once.
# tissue_integer <- as.integer(anno_col[names(clustersBestK),]$Tissue)
# tissue_integer[which(tissue_integer==11)] <- 4
# tissue_integer[which(tissue_integer==12)] <- 7
# writeMat("../data/tissue.mat", tissue = tissue_integer)

annotations <- data.frame(Cluster = as.factor(clustersBestK),
                          Tissue = anno_col[names(clustersBestK),]$Tissue,
                          COCA = anno_col[names(clustersBestK),]$COCA)

our_cl_palette <-  c(
  "1" = "#B7BF10" , "2" = "#af95a6", "3" =  "#D50032", "4" = "#B7BF10",
  "5" = "#E47113", "6" = "#64A70B", "7" = "#f7dea3", "8" = "#000000",
  "9" = "#0072ce"
)
legend_our_clusters <- Legend(title = "Clusters",
                               at = c("1", "2", "3", "4", "5",
                                      "6", "7", "8", "9"),
                               legend_gp = gpar(fill=our_cl_palette),
                               labels_gp = gpar(fontsize = 12),
                               title_gp = gpar(fontsize = 12),
                               direction = "vertical",
                               nrow = 2,
                               grid_height = unit(5, "mm"))
coca_palette <- c("#E89CAE", "#6CACE4", "#B7BF10", "#F1BE48", "#af95a6",
                  "#0072ce", "#E47113", "#93328E", "#D50032", "#85b09A",
                          "#00B0B9","#003C71")

names(coca_palette) <- c("1 - LUAD-enriched", "2 - Squamous-like",
                         "3 - BRCA/Luminal", "4 - BRCA/Basal", "5 - KIRC",
                         "6 - UCEC", "7 - COAD/READ", "8 - BLCA", "9 - OV",
                         "10 - GBM", "11 - small-various", "12 - small-various")

legend_coca_hoadley <- Legend(title = "COCA (Hoadley et al., 2014)",
                              at = c("1 - LUAD-enriched", "2 - Squamous-like",
                                    "3 - BRCA/Luminal", "4 - BRCA/Basal",
                                    "5 - KIRC", "6 - UCEC", "7 - COAD/READ",
                                    "8 - BLCA", "9 - OV", "10 - GBM",
                                    "11 - small-various",
                                    "12 - small-various"),
                              legend_gp = gpar(fill=coca_palette),
                              labels_gp = gpar(fontsize = 12),
                              title_gp = gpar(fontsize = 12),
                              direction = "vertical",
                              nrow = 7,
                              grid_height = unit(5, "mm"))

tissue_palette <- c("#93328E", "#B7BF10", "#4E5B31", "#6CACE4", "#af95a6",
                    "#E89CAE", "#8A1538", "#D50032", "#3F2A56", "#0072ce")
names(tissue_palette) <- c("BLCA", "BRCA", "COAD", "HNSC", "KIRC", "LUAD",
                           "LUSC", "OV", "READ", "UCEC")
legend_tissue <- Legend(title = "Tissue",
                        at = as.character(unique(annotations$Tissue)),
                        legend_gp = gpar(fill=tissue_palette),
                        labels_gp = gpar(fontsize = 12),
                        title_gp = gpar(fontsize = 12),
                        direction = "vertical",
                        grid_height = unit(5, "mm"),
                        nrow = 2
)

Hclusters <-
  rowAnnotation(
    df = annotations[sortedClusters$ix,],
    col = list("COCA" = coca_palette,
               "Cluster" = our_cl_palette,
               "Tissue" = tissue_palette),
    name = "Final clusters",
    show_legend = FALSE,
    show_annotation_name = TRUE)

legend_dataset <- Legend(
  title = "Weighted similarity",
  col_fun = myBlues,
  labels_gp = gpar(fontsize = 12),
  title_gp = gpar(fontsize = 12),
  direction = "horizontal",
  # grid_height = unit(5, "mm"),
)

HwMatrix <-
  Heatmap(
    KM[sortedClusters$ix, sortedClusters$ix, bestK - 1],
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    name = "Weighted kernel",
    col = myBlues,
    # row_split = sortedClusters$x,
    # column_split = sortedClusters$x,
    right_annotation = Hclusters,
    show_heatmap_legend = FALSE,
    width = unit(10, "cm"),
    height = unit(10, "cm"))

packed_legends <- packLegend(
  legend_dataset,
  legend_our_clusters,
  legend_tissue,
  legend_coca_hoadley,
  max_width= unit(10, "cm"),
  direction = "vertical",
  column_gap = unit(1, "cm"),
  row_gap = unit(1, "cm"))


pdf(
  paste0("../figures/unsupervised_weightedKernel_", bestK, "clusters", var_sel,
         ".pdf"),
  height = 5,
  width = 10
)
draw(
  HwMatrix,
  annotation_legend_list = packed_legends,
  annotation_legend_side = "right"
)
dev.off()

save(clustersBestK,
     file = paste0("../unsupervised-integration-output/", bestK,
                   "combined_clusters_", var_sel,".RData"))

############################## Plot weight matrix ##############################

legend_our_clusters <- Legend(title = "Clusters",
                              at = c("1", "2", "3", "4", "5",
                                     "6", "7", "8", "9", "10", "11",
                                     "12", "13", "14", "15"),
                              legend_gp = gpar(fill=our_cl_palette),
                              labels_gp = gpar(fontsize = 12),
                              title_gp = gpar(fontsize = 12),
                              direction = "vertical",
                              nrow = 2,
                              grid_height = unit(5, "mm"))

weight_palette<- colorRamp2(c(0,1), colors = c("white", "#0072ce"))

weightsBestK <- theta[bestK-1,,]

legend_weights <- Legend(
  col_fun = weight_palette,
  title = "Weight",
  labels_gp = gpar(fontsize = 12),
  title_gp = gpar(fontsize = 12),
  direction = "horizontal",
)

Hweights <- Heatmap(
  weightsBestK[sortedClusters$ix,],
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  col = weight_palette,
  right_annotation = Hclusters,
  show_heatmap_legend = FALSE,
  heatmap_width = unit(4, "cm")
)

packed_legends <- packLegend(
  legend_weights,
  legend_coca_hoadley,
  legend_our_clusters,
  legend_tissue,
  max_height= unit(15, "cm"),
  direction = "vertical",
  column_gap = unit(1, "cm"),
  row_gap = unit(1, "cm"))

pdf(paste0("../figures/unsupervised_weights_", bestK, var_sel,".pdf"),
    height = 5,
    width = 6)
draw(
  Hweights,
  annotation_legend_list = packed_legends,
  annotation_legend_side = "right"
)
dev.off()

############################## Coincidence matrix ##############################

### KLIC vs COCA (Hoadley et al.)

klicK <- bestK
cocaK <- 12

coincidences <- matrix(NA, cocaK, klicK)

rownames(coincidences) <- paste("COCA cluster", 1:cocaK)
colnames(coincidences) <- paste("Cluster", 1:klicK)

COCALabels <- as.numeric(substr(annotations$COCA, start = 1, stop=2))
names(COCALabels) <- rownames(annotations)

KLICLabels <- annotations$Cluster
names(KLICLabels) <- rownames(annotations) 

for(i in 1:cocaK){
  whichCOCA_i <- names(COCALabels)[which(COCALabels == i)]
  cat("Number of elements of COCA in cluster",
      i, ":", length(whichCOCA_i), '\n')
  for(j in 1:klicK){
    whichKLIC_j <- names(KLICLabels)[which(KLICLabels == j)]
    cat("Number of elements of KLIC in cluster",
        j, ":", length(whichKLIC_j), "\n")
    coincidences[i,j] <- sum(whichCOCA_i %in% whichKLIC_j)
    cat("Intersection: ", coincidences[i,j], "\n")
  }
}

col_fun = colorRamp2(c(0, max(coincidences)), c("white", "#D50032"))

pdf(paste0("../figures/unsupervised_klic_coca_coincidences_", klicK,"_", cocaK,
           var_sel, ".pdf"),
    height = 6, 
    width = 9)
Heatmap(
  coincidences,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%d", as.integer(coincidences[i, j])), x, y,
              gp = gpar(fontsize = 16))
  },
  col = col_fun,
  row_names_gp = gpar(fontsize = 14),
  column_names_gp = gpar(fontsize = 14),
  heatmap_legend_param = list(
    title = "Coincidences",
    title_gp = gpar(fontsize = 16),
    labels_gp = gpar(fontsize = 16),
    direction = "vertical",
    nrow = 1,
    title_position = "lefttop-rot",
    legend_height = unit(6, "cm")
  ),
  width = unit(bestK, "cm")
)
dev.off()

 ### KLIC vs Tissue

klicK <- bestK
tissueK <- 10 # Note that we have removed two tissues 

coincidences <- matrix(NA, tissueK, klicK)

rownames(coincidences) <- names(which(table(annotations$Tissue) != 0))
colnames(coincidences) <- paste("Cluster", 1:klicK)

TissueLabels <- annotations$Tissue
names(TissueLabels) <- rownames(annotations)

KLICLabels <- annotations$Cluster
names(KLICLabels) <- rownames(annotations) 

for(i in rownames(coincidences)){
  whichTissue_i <- names(TissueLabels)[which(TissueLabels == i)]
  cat("Number of samples of tissue", i, ":", length(whichTissue_i), '\n')
  for(j in 1:klicK){
    whichKLIC_j <- names(KLICLabels)[which(KLICLabels == j)]
    cat("Number of elements of KLIC in cluster", j, ":", length(whichKLIC_j),
        "\n")
    coincidences[i,j] <- sum(whichTissue_i %in% whichKLIC_j)
    cat("Intersection: ", coincidences[i,j], "\n")
  }
}

col_fun = colorRamp2(c(0, max(coincidences)), c("white", "#D50032"))

pdf(paste0("../figures/unsupervised_klic_tissue_coincidences_",
           klicK, "_", tissueK, var_sel, ".pdf"),
    height = 6, 
    width = 8)
Heatmap(
  coincidences,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%d", as.integer(coincidences[i, j])), x, y,
              gp = gpar(fontsize = 16))
  },
  col = col_fun,
  row_names_gp = gpar(fontsize = 14),
  column_names_gp = gpar(fontsize = 14),
  heatmap_legend_param = list(
    title = "Coincidences",
    title_gp = gpar(fontsize = 16),
    labels_gp = gpar(fontsize = 16),
    direction = "vertical",
    nrow = 1,
    title_position = "lefttop-rot",
    legend_height = unit(6, "cm")
  ),
  width = unit(bestK, "cm")
)
dev.off()

################################ Plot silhouette ###############################

pdf(paste0("../figures/unsupervised_pancan_silhouette", var_sel,".pdf"),
    height = 3.75,
    width = 6
    )
plot(2:maxK,
     maxSil$silhouette,
     type = "b",
     xlab = "Number of clusters",
     ylab = "Silhouette")
dev.off()
