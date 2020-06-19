##################### Pancancer study (Hoadley et al. 2014) ####################
############################ PSM of copy number data ###########################

rm(list=ls())

library(ComplexHeatmap)
library(klic)
library(mclust)

################################## Compute PSM #################################

mcmc_samples <- read.csv("../mdi/methylation_50max.csv")
# mcmc_samples <- mcmc_samples[seq.int(from = 1001, to = 2000), 2:2422]
# colnames(mcmc_samples) <- sub("Dataset1_", "", colnames(mcmc_samples))
# 
# psm_methylation <- matrix(0, dim(mcmc_samples)[2], dim(mcmc_samples)[2])
# n_clusters <- rep(0, dim(mcmc_samples)[1])
# for(i in 1:dim(mcmc_samples)[1]){
#   cat(i, "\n")
#   sample_i <-mcmc_samples[i,]
#   names(sample_i) <- NULL
#   sample_i <- unlist(sample_i)
#   n_clusters[i] <- length(unique(sample_i))
#   for(j in unique(sample_i)){
#     psm_methylation <- psm_methylation + crossprod(mcmc_samples[i,] == j)
#   }
# }
# psm_methylation <- psm_methylation / dim(mcmc_samples)[1]
# coph_corr_methylation <- copheneticCorrelation(psm_methylation)

################################## Cluster PSM #################################

load("../data/samples.RData")

# table(n_clusters)
# 
# parameters <- list(
#   cluster_count =
#     as.numeric(names(table(n_clusters))[which.max(table(n_clusters))]))
# psm_methylation <-
#   spectrumShift(psm_methylation, verbose = TRUE, shift = 0.000001)
# kkm_output <- kkmeans(psm_methylation, parameters)
# clusters <- kkm_output$clustering
# names(clusters) <- gsub("\\.", "-", rownames(psm_methylation))
# 
# ari <- adjustedRandIndex(clusters, anno_col[names(clusters),]$Tissue) 
# save(psm_methylation, coph_corr_methylation, n_clusters, clusters, ari,
#      file = "../mdi/psm_methylation.RData")
load("../mdi/psm_methylation.RData")

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

rownames(psm_methylation) <- rownames(psm_methylation) <-
  gsub("\\.", "-", rownames(psm_methylation))
annotations_methylation <- anno_col[rownames(anno_col) %in%
                                      rownames(psm_methylation),]
annotations_methylation <- annotations_methylation[rownames(psm_methylation),]

Hclusters <-
  rowAnnotation(Tissue = as.character(annotations_methylation$Tissue),
                name = "Final clusters",
                annotation_legend_param =
                my_anno_legend_param("Tissue",nrow=2),
                show_annotation_name = FALSE,
                col = list(Tissue =  palette12))

Hpsm <- Heatmap(psm_methylation,
                col = myBlues,
                show_column_names = FALSE,
                show_row_names = FALSE,
                heatmap_legend_param = my_heatmap_legend_param,
                name = "Similarity",
                right_annotation = Hclusters,
                show_column_dend = FALSE,
                show_row_dend = FALSE
)

png("../figures/psm_mdi_methylation.png",
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

mass_parameter <- mcmc_samples[,1]

png("../figures/mass_parameter_methylation_chain.png",
height = 400, width = 400)
plot(seq.int(from=5001, to=10000, by = 5),
     mass_parameter[1001:2000],
     type = 'l', xlab = "Iteration", ylab = "Mass parameter")
dev.off()

png("../figures/mass_parameter_methylation_posterior.png",
    height = 400, width = 400)
hist(mass_parameter[1001:2000],
     breaks = 50,
     main = "",
     xlab = "Mass parameter")
dev.off()

png("../figures/n_clusters_methylation_chain.png",
    height = 400, width = 400)
plot(seq.int(from=5001, to=10000, by = 5),
     n_clusters,
     type = 'l',
     xlab = "Iteration",
     ylab = "Number of clusters")
dev.off()

png("../figures/n_clusters_methylation_posterior.png",
    height = 400, width = 400)
hist(n_clusters,
     # breaks = 50,
     main = "",
     xlab = "Number of clusters")
dev.off()
