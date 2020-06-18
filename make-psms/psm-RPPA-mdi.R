##################### Pancancer study (Hoadley et al. 2014) ####################
############################ PSM of copy number data ###########################

rm(list=ls())

library(ComplexHeatmap)
library(klic)
library(mclust)

################################# Compute PSM ##################################

# samples <- read.csv("../mdi/RPPA_50max.csv")
# samples <- samples[seq.int(from = 1001, to = 2000), 2:2422]
# colnames(samples) <- sub("Dataset1_", "", colnames(samples))
# 
# psm_rppa <- matrix(0, dim(samples)[2], dim(samples)[2])
# n_clusters <- rep(0, dim(samples)[1])
# for(i in 1:dim(samples)[1]){
#   cat(i, "\n")
#   sample_i <-samples[i,]
#   names(sample_i) <- NULL
#   sample_i <- unlist(sample_i)
#   n_clusters[i] <- length(unique(sample_i))
#   for(j in unique(sample_i)){
#     psm_rppa <- psm_rppa + crossprod(samples[i,] == j)
#   }
# }
# psm_rppa <- psm_rppa / dim(samples)[1]
# coph_corr_rppa <- copheneticCorrelation(psm_rppa)

################################## Cluster PSM #################################

load("../data/samples.RData")

# table(n_clusters)
# 
# parameters <- list(
#   cluster_count =
#     as.numeric(names(table(n_clusters))[which.max(table(n_clusters))]))
# psm_rppa <- spectrumShift(psm_rppa, verbose = TRUE, shift = 0.000001)
# kkm_output <- kkmeans(psm_rppa, parameters)
# clusters <- kkm_output$clustering
# names(clusters) <- gsub("\\.", "-", rownames(psm_rppa))
# 
# ari <- adjustedRandIndex(clusters, anno_col[names(clusters),]$Tissue) 
# save(psm_rppa, coph_corr_rppa, n_clusters, clusters, ari,
#      file = "../mdi/psm_rppa.RData")
load("../mdi/psm_rppa.RData")

################################### Plot PSM ###################################

palette25 = c("#B7BF10", # Light green
              "#4E5B31", # Core green
              "#115E67", # Core Cambridge blue
              "#85b09A", # Light Cambridge blue
              "#0072ce", # Core blue
              "#6CACE4", # Light blue
              "#E89CAE", # Light pink
              "#af95a6", # Light purple
              "#8C77A3", # Modified core purple
              "#D50032", # Core red
              "#E87722", # Core orange
              "#F1BE48", # Light yellow
              "#edf373", # Light green
              "#adbe86", # Core green
              "#57d4e3", # Core Cambridge blue
              "#c2d7cc", # Light Cambridge blue
              "#66bbff", # Core blue
              "#b5d5f1", # Light blue
              "#f3cdd6", # Light pink
              "#d7cad2", # Light purple
              "#c5bbd1", # Modified core purple
              "#ff6a8d", # Core red
              "#f3bb90", # Core orange
              "#f7dea3", # Light yellow
              "#000000") # Black

palette12 <- palette25[sample(1:12,12)]
names(palette12) <- names(table(anno_col$Tissue))

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

rownames(psm_rppa) <- rownames(psm_rppa) <- gsub("\\.", "-", rownames(psm_rppa))
annotations_rppa <- anno_col[rownames(anno_col) %in% rownames(psm_rppa),]
annotations_rppa <- annotations_rppa[rownames(psm_rppa),]

Hclusters <-
  rowAnnotation(Tissue = as.character(annotations_rppa$Tissue),
                name = "Final clusters",
                annotation_legend_param = my_anno_legend_param("Tissue",nrow=2),
                show_annotation_name = FALSE,
                col = list(Tissue =  palette12))

Hpsm <- Heatmap(psm_rppa,
                col = myBlues,
                show_column_names = FALSE,
                show_row_names = FALSE,
                heatmap_legend_param = my_heatmap_legend_param,
                name = "Similarity",
                right_annotation = Hclusters,
                show_column_dend = FALSE,
                show_row_dend = FALSE
                )

png("../figures/psm_mdi_rppa.png",
  height = 572,
  width = 520,
)
draw(
  Hpsm,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)
dev.off()

save(palette12, file = "../data/tumour_colours.RData")
