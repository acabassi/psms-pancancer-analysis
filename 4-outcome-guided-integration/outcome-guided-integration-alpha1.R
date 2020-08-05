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

var_sel <- "_alpha1"
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

################################# Load weights #################################

if(var_sel == ""){
  var_sel = "_noVarSel"
}
weights <-
  as.vector(
    readMat(
      paste0("../outcome-guided-integration-output/outcome_guided_integration",
             var_sel, "_weights.mat"))$beta)

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
weighted_kernel <- array(NA, c(N, N))
dimnames(weighted_kernel) <- list(rownames(psm_cn), rownames(psm_cn))

weighted_kernel <- weights[1]*psm_cn + weights[2]*psm_mrna +
  weights[3]*psm_methylation + weights[4]*psm_mirna + weights[5]*psm_rppa

################################# Integration ##################################
# Use kernel k-means to integrate the datasets

# Initialise array of kernel matrices
maxK <- 50
clLabels <- array(NA, c(maxK - 1, N))

### Each of these must is run separately on the HPC ###
parameters <- list()
parameters$iteration_count <-
  500 # set the maximum number of iterations

for(i in 2:50){
  # Use kernel k-means with K=i to find weights and cluster labels
  parameters$cluster_count <- i # set the number of clusters K
  kkm <- kkmeans(weighted_kernel, parameters, seed = 1)
  
  # Save cluster labels
  clLabels[i - 1, ] <- kkm$clustering
}

############################### Maximise silhouette ############################

## This is silly but I need it as input of the maximiseSilhouette() function
repeated_weighted_kernel <- array(NA, c(N, N, maxK-1))
for(i in 1:(maxK-1)){
  repeated_weighted_kernel[,,i] <- weighted_kernel 
}

maxSil <- maximiseSilhouette(repeated_weighted_kernel, clLabels, maxK = maxK)

save(
  weighted_kernel,
  maxSil,
  file = paste0(
    "../outcome-guided-integration-output/outcome_guided_output_all",
    var_sel,".RData")
)


bestK <- maxSil$K
inds <- rownames(weighted_kernel)
clustersBestK <- clLabels[bestK - 1, ]
names(clustersBestK) <- rownames(weighted_kernel)

################################## Plot output #################################

# Let's just plot the clusters:
table(clustersBestK)
sortedClusters <- sort(clustersBestK, index.return = T)

myBlues <- colorRamp2(c(0,1), colors = c("white", "#003C71"))
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

annotations <- data.frame(Cluster = as.factor(clustersBestK),
                          Tissue = anno_col[names(clustersBestK),]$Tissue,
                          COCA = anno_col[names(clustersBestK),]$COCA)
our_cl_palette <-  c(
  "1" = "#af95a6" , "2" = "#af95a6", "3" =  "#F1BE48", "4" = "#0072ce",
  "5" = "#0072ce", "6" = "#93328E", "7" = "#B7BF10", "8" = "#B7BF10",
  "9" = "#B7BF10", "10" ="#af95a6", "11" = "#B7BF10", "12" = "#6CACE4",
  "13" = "#6CACE4", "14" = "#E89CAE", "15" = "#0072ce", "16" = "#af95a6",
  "17" = "#B7BF10", "18" = "#D50032", "19" = "#B7BF10",
  "20" = "#B7BF10", "21" = "#8A1538", "22" = "#E47113"
)
legend_our_clusters <- Legend(title = "Clusters",
                              at = c("1", "2", "3", "4", "5",
                                     "6", "7", "8", "9", "10", "11",
                                     "12", "13", "14", "15", "16", "17", "18",
                                     "19", "20", "21", "22"),
                              legend_gp = gpar(fill=our_cl_palette),
                              labels_gp = gpar(fontsize = 12),
                              title_gp = gpar(fontsize = 12),
                              direction = "vertical",
                              nrow = 4,
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
                              nrow = 6,
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
    weighted_kernel[sortedClusters$ix, sortedClusters$ix],
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    name = "Weighted kernel",
    col = myBlues,
    # row_split = sortedClusters$x,
    # column_split = sortedClusters$x,
    right_annotation = Hclusters,
    show_heatmap_legend = FALSE,
    width = unit(13, "cm"),
    height = unit(13, "cm"))

packed_legends <- packLegend(
  legend_dataset,
  legend_our_clusters,
  legend_tissue,
  legend_coca_hoadley,
  max_height = unit(13,"cm"),
  max_width= unit(10, "cm"),
  direction = "vertical",
  column_gap = unit(1, "cm"),
  row_gap = unit(1, "cm"))


pdf(
  paste0("../figures/outcomeguided_weightedKernel_", bestK, "clusters",
         var_sel, ".pdf"),
  height = 6,
  width = 10
)
draw(
  HwMatrix,
  annotation_legend_list = packed_legends,
  annotation_legend_side = "right"
)
dev.off()

save(clustersBestK,
     file = paste0("../outcome-guided-integration-output/", bestK,
                   "_combined_clusters_", var_sel,".RData"))

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

pdf(paste0("../figures/outcomeguided_klic_coca_coincidences_", klicK,"_", cocaK,
            var_sel, ".pdf"),
    height = 6, 
    width = 16)
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

pdf(paste0("../figures/outcomeguided_klic_tissue_coincidences_",
           klicK, "_", tissueK, var_sel, ".pdf"),
    height = 6, 
    width = 14)
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

pdf(paste0("../figures/outcomeguided_pancan_silhouette_", var_sel,".pdf"),
    height = 3.75,
    width = 6
    )
plot(2:maxK,
     maxSil$silhouette,
     type = "b",
     xlab = "Number of clusters",
     ylab = "Silhouette")
dev.off()

