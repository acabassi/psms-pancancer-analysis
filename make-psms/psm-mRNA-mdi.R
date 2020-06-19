##################### Pancancer study (Hoadley et al. 2014) ####################
############################ PSM of copy number data ###########################

rm(list=ls())

library(ComplexHeatmap)
library(klic)
library(mclust)

#### Check that maximum number of clusters has been set correctly in mdipp #####

# samples_max100 <- read.csv("../mdi/mRNA_100max.csv")
# n_clusters_max100 <- ari<- rep(0, 1000)
# for(i in 1:dim(samples_max100)[1]){
#   sample_i <-samples_max100[i,]
#   names(sample_i) <- NULL
#   sample_i <- unlist(sample_i)
#   n_clusters_max100[i] <- length(unique(sample_i))
#   ari[i]<-adjustedRandIndex(sample_i, samples_max100[1000,])
# }
# plot(1:1000, n_clusters_max100)
# plot(ari)

################################# Compute PSM ##################################

mcmc_samples <- read.csv("../mdi/mRNA_200max.csv")
# samples <- samples[seq.int(from = 1001, to = 2000), 2:2422]
# colnames(samples) <- sub("Dataset1_", "", colnames(samples))
# 
# psm_mrna <- matrix(0, dim(samples)[2], dim(samples)[2])
# n_clusters <- rep(0, dim(samples)[1])
# for(i in 1:dim(samples)[1]){
#   cat(i, "\n")
#   sample_i <-samples[i,]
#   names(sample_i) <- NULL
#   sample_i <- unlist(sample_i)
#   n_clusters[i] <- length(unique(sample_i))
#   for(j in unique(sample_i)){
#     psm_mrna <- psm_mrna + crossprod(samples[i,] == j)
#   }
# }
# psm_mrna <- psm_mrna / dim(samples)[1]
# coph_corr_mrna <- copheneticCorrelation(psm_mrna)

################################## Cluster PSM #################################

load("../data/samples.RData")

# table(n_clusters)
# 
# parameters <- list(
#   cluster_count =
#     as.numeric(names(table(n_clusters))[which.max(table(n_clusters))]))
# psm_mrna <- spectrumShift(psm_mrna, verbose = TRUE, shift = 0.000001)
# kkm_output <- kkmeans(psm_mrna, parameters)
# clusters <- kkm_output$clustering
# names(clusters) <- gsub("\\.", "-", rownames(psm_mrna))
# 
# ari <- adjustedRandIndex(clusters, anno_col[names(clusters),]$Tissue) 
# save(psm_mrna, coph_corr_mrna, n_clusters, clusters, ari,
#      file = "../mdi/psm_mrna.RData")
load("../mdi/psm_mrna.RData")

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

rownames(psm_mrna) <- rownames(psm_mrna) <-
  gsub("\\.", "-", rownames(psm_mrna))
annotations_mrna <- anno_col[rownames(anno_col) %in% rownames(psm_mrna),]
annotations_mrna <- annotations_mrna[rownames(psm_mrna),]

Hclusters <-
  rowAnnotation(Tissue = as.character(annotations_mrna$Tissue),
                name = "Final clusters",
                annotation_legend_param = my_anno_legend_param("Tissue",nrow=2),
                show_annotation_name = FALSE,
                col = list(Tissue =  palette12))

Hpsm <- Heatmap(psm_mrna,
                col = myBlues,
                show_column_names = FALSE,
                show_row_names = FALSE,
                heatmap_legend_param = my_heatmap_legend_param,
                name = "Similarity",
                right_annotation = Hclusters,
                show_column_dend = FALSE,
                show_row_dend = FALSE
)

png("../figures/psm_mdi_mrna.png",
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

png("../figures/mass_parameter_mrna_chain.png",
    height = 400, width = 400)
plot(seq.int(from=5001, to=10000, by = 5),
     mass_parameter[1001:2000],
     type = 'l', xlab = "Iteration", ylab = "Mass parameter")
dev.off()

png("../figures/mass_parameter_mrna_posterior.png",
    height = 400, width = 400)
hist(mass_parameter[1001:2000],
     breaks = 50,
     main = "",
     xlab = "Mass parameter")
dev.off()

png("../figures/n_clusters_mrna_chain.png",
    height = 400, width = 400)
plot(seq.int(from=5001, to=10000, by = 5),
     n_clusters,
     type = 'l',
     xlab = "Iteration",
     ylab = "Number of clusters")
dev.off()

png("../figures/n_clusters_mrna_posterior.png",
    height = 400, width = 400)
hist(n_clusters,
     # breaks = 50,
     main = "",
     xlab = "Number of clusters")
dev.off()
