##################### Pancancer study (Hoadley et al. 2014) ####################
############################ PSM of copy number data ###########################

rm(list=ls())

library(ComplexHeatmap)
library(klic)
library(mclust)

#### Check that maximum number of clusters has been set correctly in mdipp #####

# samples_max200 <- read.csv("../mdi/CN_200max.csv")
# n_clusters_max200 <- rep(0, 1500)
# for(i in 1:dim(samples_max200)[1]){
#   sample_i <-samples_max200[i,]
#   names(sample_i) <- NULL
#   sample_i <- unlist(sample_i)
#   n_clusters_max200[i] <- length(unique(sample_i))
# }
# plot(1:1500, n_clusters_max200)

################################## Compute PSM #################################

# samples <- read.csv("../mdi/CN.csv")
# samples <- samples[seq.int(from = 1001, to = 2000), 2:2422]
# colnames(samples) <- sub("Dataset1_", "", colnames(samples))
# 
# psm_cn <- matrix(0, dim(samples)[2], dim(samples)[2])
# n_clusters <- rep(0, dim(samples)[1])
# for(i in 1:dim(samples)[1]){
#   cat(i, "\n")
#   sample_i <-samples[i,]
#   names(sample_i) <- NULL
#   sample_i <- unlist(sample_i)
#   n_clusters[i] <- length(unique(sample_i))
#   for(j in unique(sample_i)){
#     psm_cn <- psm_cn + crossprod(samples[i,] == j)
#   }
# }
# psm_cn <- psm_cn / dim(samples)[1]
# coph_corr_cn <- copheneticCorrelation(psm_cn)

################################## Cluster PSM #################################

load("../data/samples.RData")


# table(n_clusters)
# 
# parameters <- list(
#   cluster_count =
#     as.numeric(names(table(n_clusters))[which.max(table(n_clusters))]))
# psm_cn <- spectrumShift(psm_cn, verbose = TRUE, shift = 0.000001)
# kkm_output <- kkmeans(psm_cn, parameters)
# clusters <- kkm_output$clustering
# names(clusters) <- gsub("\\.", "-", rownames(psm_cn))
# 
# ari <- adjustedRandIndex(clusters, anno_col[names(clusters),]$Tissue) 
# save(psm_cn, coph_corr_cn, n_clusters, clusters, ari,
#      file = "../mdi/psm_cn.RData")
load("../mdi/psm_cn.RData")

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

png("../figures/psm_mdi_cn.png",
    height = 572,
    width = 520
)
draw(
  Hpsm,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)
dev.off()


