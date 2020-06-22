##################### Pancancer study (Hoadley et al. 2014) ####################
############################ PSM of copy number data ###########################

rm(list=ls())

library(ComplexHeatmap)
library(klic)
library(mclust)

################################# Select chain #################################

chain <- 1 # Must be an integer between 1 and 5

#### Check that maximum number of clusters has been set correctly in mdipp #####

mcmc_samples <- read.csv(paste0("../mdi/CN_250max_chain",chain,".csv"))
n_clusters <- rep(0, dim(mcmc_samples)[1])
for(i in 1:dim(mcmc_samples)[1]){
  sample_i <- mcmc_samples[i,]
  names(sample_i) <- NULL
  sample_i <- unlist(sample_i)
  n_clusters[i] <- length(unique(sample_i))
}
plot(n_clusters)

################################## Compute PSM #################################

mcmc_samples <- mcmc_samples[1001:2000, 2:2422]
colnames(mcmc_samples) <- sub("Dataset1_", "", colnames(mcmc_samples))
Rcpp::sourceCpp("makePSM.cpp")
psm_cn <- makePSM(as.matrix(mcmc_samples))
rownames(psm_cn) <- colnames(psm_cn) <- colnames(mcmc_samples)
coph_corr_cn <- copheneticCorrelation(psm_cn)

################################## Cluster PSM #################################

load("../data/samples.RData")

table(n_clusters)

parameters <- list(
  cluster_count =
    as.numeric(names(table(n_clusters))[which.max(table(n_clusters))]))
psm_cn <- spectrumShift(psm_cn, verbose = TRUE, shift = 0.000001)
kkm_output <- kkmeans(psm_cn, parameters)
clusters <- kkm_output$clustering
names(clusters) <- gsub("\\.", "-", rownames(psm_cn))

ari <- adjustedRandIndex(clusters, anno_col[names(clusters),]$Tissue)
save(psm_cn, coph_corr_cn, n_clusters, clusters, ari,
     file = paste0("../mdi/psm_cn_chain", chain,".RData"))
load(paste0("../mdi/psm_cn_chain", chain,".RData"))

################################### Plot PSM ###################################

load("../data/tumour_colours.RData")

myBlues <-
  colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)

my_anno_legend_param <- function(name, nrow=1) {
  return(list(
    title = name,
    labels_gp = gpar(fontsize = 18),
    title_gp = gpar(fontsize = 22),
    nrow = nrow,
    grid_height = unit(5, "mm"),
    grid_width = unit(5, "mm")
  ))
}

my_heatmap_legend_param <- list(labels_gp = gpar(fontsize = 18),
                                title_gp = gpar(fontsize = 22),
                                direction = "horizontal",
                                nrow = 1,
                                at = c(0,0.5,1),
                                legend_width = unit(4, "cm"))

rownames(psm_cn) <- rownames(psm_cn) <- gsub("\\.", "-", rownames(psm_cn))
annotations_cn <- anno_col[rownames(anno_col) %in% rownames(psm_cn),]
annotations_cn <- annotations_cn[rownames(psm_cn),]

Hclusters <-
  rowAnnotation(Tissue = as.character(annotations_cn$Tissue),
                name = "Final clusters",
                annotation_legend_param =
                  my_anno_legend_param("Tissue",nrow=2),
                show_annotation_name = FALSE,
                col = list(Tissue =  palette12))

Hpsm <- Heatmap(psm_cn,
                col = myBlues,
                show_column_names = FALSE,
                show_row_names = FALSE,
                heatmap_legend_param = my_heatmap_legend_param,
                name = "Similarity",
                right_annotation = Hclusters,
                show_column_dend = FALSE,
                show_row_dend = FALSE
)

png(paste0("../figures/psm_mdi_cn_chain", chain,".png"),
    height = 572,
    width = 520
)
draw(
  Hpsm,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)
dev.off()

########################### Convergence  assessment ############################

mcmc_samples <- read.csv(paste0("../mdi/CN_250max_chain",chain,".csv"))
mass_parameter <- mcmc_samples[,1]

png(paste0("../figures/mass_parameter_cn_chain", chain,".png"),
    height = 400, width = 400)
plot(mass_parameter,
     type = 'l', xlab = "Iteration", ylab = "Mass parameter")
dev.off()

png(paste0("../figures/mass_parameter_cn_posterior_chain", chain, ".png"),
    height = 400, width = 400)
hist(mass_parameter[1001:2000],
     breaks = 50,
     main = "",
     xlab = "Mass parameter")
dev.off()

png(paste0("../figures/n_clusters_cn_chain", chain,".png"),
    height = 400, width = 400)
plot(n_clusters,
     type = 'l',
     xlab = "Iteration",
     ylab = "Number of clusters")
dev.off()

png(paste0("../figures/n_clusters_cn_posterior", chain,".png"),
    height = 400, width = 400)
hist(n_clusters[1001:2000],
     breaks = 6,
     main = "",
     xlab = "Number of clusters")
dev.off()
 
