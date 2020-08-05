### Combine PSMs for different MCMC chains

rm(list=ls())

library(ComplexHeatmap)
library(klic)
library(R.matlab)

n_chains <- 5
noCuda <- "_noCuda"
n_iterations <- 50000

# for(layer in c("CN", "mRNA", "methylation", "miRNA", "RPPA")){ 
layer <- "RPPA"
# for(var_sel in c("_alpha1", "_alpha05")){ # , "_alpha01", "_alpha05", "_alpha1"
var_sel <- "_alpha1"

if(layer %in% c("CN", "mRNA")){
  max_clusters = 250
}else if(layer %in% c("RPPA", "methylation")){
  max_clusters = 50
}else if(layer == "miRNA"){
  max_clusters = 100
}

if(noCuda == "_noCuda" & layer == "mRNA"){
  noCuda <- "_temp"
}

load(paste0(
  "../mdi/psms-rdata/psm_", layer, "_chain", 1, var_sel, noCuda, ".RData"))
psm_cn1 <- psm_cn; rm(psm_cn)
load(paste0(
  "../mdi/psms-rdata/psm_", layer, "_chain", 2, var_sel, noCuda, ".RData"))
psm_cn2 <- psm_cn; rm(psm_cn)
load(paste0(
  "../mdi/psms-rdata/psm_", layer, "_chain", 3, var_sel, noCuda, ".RData"))
psm_cn3 <- psm_cn; rm(psm_cn)
load(paste0(
  "../mdi/psms-rdata/psm_", layer, "_chain", 4, var_sel, noCuda, ".RData"))
psm_cn4 <- psm_cn; rm(psm_cn)
load(paste0(
  "../mdi/psms-rdata/psm_", layer, "_chain", 5, var_sel, noCuda, ".RData"))
psm_cn5 <- psm_cn; rm(psm_cn)

psm_list <- list(psm_cn1, psm_cn2, psm_cn3, psm_cn4, psm_cn5)
average_psm <- Reduce("+", psm_list) / length(psm_list)

file_name <- paste0("psm_", layer, "_average_", var_sel, noCuda)
save(average_psm, file = paste0("../mdi/psms-rdata/", file_name, ".RData"))
writeMat(paste0("../mdi/psms-mat/", file_name, ".mat"),
         average_psm = average_psm)
writeMat(paste0("../mdi/psms-mat/", file_name, "_subset.mat"),
         average_psm = average_psm[1:800,1:800])

############################### Plot average PSM ###############################

# load("../data/samples.RData")
# load("../data/tumour_colours.RData")
# 
# myBlues <-
#   colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)
# 
# my_anno_legend_param <- function(name, nrow=1) {
#   return(list(
#     title = name,
#     labels_gp = gpar(fontsize = 18),
#     title_gp = gpar(fontsize = 22),
#     nrow = nrow,
#     grid_height = unit(5, "mm"),
#     grid_width = unit(5, "mm")
#   ))
# }
# 
# my_heatmap_legend_param <- list(labels_gp = gpar(fontsize = 18),
#                                 title_gp = gpar(fontsize = 22),
#                                 direction = "horizontal",
#                                 nrow = 1,
#                                 at = c(0,0.5,1),
#                                 legend_width = unit(4, "cm"))
# 
# rownames(average_psm) <- rownames(average_psm) <-
#   gsub("\\.", "-", rownames(average_psm))
# annotations_cn <- anno_col[rownames(anno_col) %in% rownames(average_psm),]
# annotations_cn <- annotations_cn[rownames(average_psm),]
# 
# Hclusters <-
#   rowAnnotation(Tissue = as.character(annotations_cn$Tissue),
#                 name = "Final clusters",
#                 annotation_legend_param =
#                   my_anno_legend_param("Tissue",nrow=2),
#                 show_annotation_name = FALSE,
#                 col = list(Tissue =  palette12))
# 
# Hpsm <- Heatmap(average_psm,
#                 col = myBlues,
#                 show_column_names = FALSE,
#                 show_row_names = FALSE,
#                 heatmap_legend_param = my_heatmap_legend_param,
#                 name = "Similarity",
#                 right_annotation = Hclusters,
#                 show_column_dend = FALSE,
#                 show_row_dend = FALSE
# )
# 
# pdf(paste0("../figures/psm_mdi_", layer, "_average_", var_sel, noCuda, ".pdf"),
#     height = 9,
#     width = 8
# )
# draw(
#   Hpsm,
#   heatmap_legend_side = "bottom",
#   annotation_legend_side = "bottom"
# )
# dev.off()
# 
