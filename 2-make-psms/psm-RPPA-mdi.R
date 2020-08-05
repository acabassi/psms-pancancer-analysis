##################### Pancancer study (Hoadley et al. 2014) ####################
############################ PSM of copy number data ###########################

rm(list=ls())

library(ComplexHeatmap)
library(klic)
library(mclust)

for(chain in 1:5){
  for(var_sel in c("", "_alpha01", "_alpha05", "_alpha1")){

#### Check that maximum number of clusters has been set correctly in mdipp #####
    
    mcmc_samples <- read.csv(paste0("../mdi/RPPA_50max_chain",
                                    chain,
                                    var_sel,
                                    ".csv"))
    
    n_clusters <- ari<- rep(0, )
    for(i in 1:dim(mcmc_samples)[1]){
      sample_i <-mcmc_samples[i,]
      names(sample_i) <- NULL
      sample_i <- unlist(sample_i)
      n_clusters[i] <- length(unique(sample_i))
      ari[i]<-adjustedRandIndex(sample_i, mcmc_samples[2000,])
    }
    plot(n_clusters)
    plot(ari)
    
################################# Compute PSM ##################################

    mass_parameter <- mcmc_samples[,1]
    # Remove burn-in
    mcmc_samples <- mcmc_samples[seq.int(from = 1001, to = 2000), 2:2422]
    colnames(mcmc_samples) <- sub("Dataset1_", "", colnames(mcmc_samples))
    Rcpp::sourceCpp("makePSM.cpp")
    psm_rppa <- makePSM(as.matrix(mcmc_samples))
    rownames(psm_rppa) <- colnames(psm_rppa) <- colnames(mcmc_samples)
    coph_corr_rppa <- copheneticCorrelation(psm_rppa)

################################## Cluster PSM #################################

    load("../data/samples.RData")

    table(n_clusters)

    parameters <- list(
      cluster_count =
        as.numeric(names(table(n_clusters))[which.max(table(n_clusters))]))
    psm_rppa <- spectrumShift(psm_rppa, verbose = TRUE, shift = 0.000001)
    kkm_output <- kkmeans(psm_rppa, parameters)
    clusters <- kkm_output$clustering
    names(clusters) <- gsub("\\.", "-", rownames(psm_rppa))

    ari <- adjustedRandIndex(clusters, anno_col[names(clusters),]$Tissue)
    save(psm_rppa, coph_corr_rppa, n_clusters, clusters, ari,
         file = paste0("../mdi/psm_rppa_chain", chain, var_sel, ".RData"))
    load(paste0("../mdi/psm_rppa_chain", chain, var_sel, ".RData"))

################################### Plot PSM ###################################

# To choose the tumour colours for all plots
# Do not run multiple times
#   palette25 = c("#B7BF10", # Light green
#                 "#4E5B31", # Core green
#                 "#115E67", # Core Cambridge blue
#                 "#85b09A", # Light Cambridge blue
#                 "#0072ce", # Core blue
#                 "#6CACE4", # Light blue
#                 "#E89CAE", # Light pink
#                 "#af95a6", # Light purple
#                 "#8C77A3", # Modified core purple
#                 "#D50032", # Core red
#                 "#E87722", # Core orange
#                 "#F1BE48", # Light yellow
#                 "#edf373", # Light green
#                 "#adbe86", # Core green
#                 "#57d4e3", # Core Cambridge blue
#                 "#c2d7cc", # Light Cambridge blue
#                 "#66bbff", # Core blue
#                 "#b5d5f1", # Light blue
#                 "#f3cdd6", # Light pink
#                 "#d7cad2", # Light purple
#                 "#c5bbd1", # Modified core purple
#                 "#ff6a8d", # Core red
#                 "#f3bb90", # Core orange
#                 "#f7dea3", # Light yellow
#                 "#000000") # Black
# 
# palette12 <- palette25[sample(1:12,12)]
# names(palette12) <- names(table(anno_col$Tissue))

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
    
    rownames(psm_rppa) <- rownames(psm_rppa) <- gsub("\\.", "-",
                                                     rownames(psm_rppa))
    annotations_rppa <- anno_col[rownames(anno_col) %in% rownames(psm_rppa),]
    annotations_rppa <- annotations_rppa[rownames(psm_rppa),]
    
    Hclusters <-
      rowAnnotation(Tissue = as.character(annotations_rppa$Tissue),
                    name = "Final clusters",
                    annotation_legend_param =
                      my_anno_legend_param("Tissue",nrow=2),
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
    
    pdf(paste0("../figures/psm_mdi_rppa_chain",
               chain,
               var_sel,
               ".pdf"),
      height = 572,
      width = 520,
    )
    draw(
      Hpsm,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom"
    )
    dev.off()
    
    # save(palette12,
    #      file = paste0("../data/tumour_colours_chain",
    #                    chain,
    #                    var_sel,
    #                    ".RData"))
    
########################### Convergence  assessment ############################
    
    pdf(paste0("../figures/mass_parameter_rppa_chain",
               chain,
               var_sel,
               ".pdf"),
        height = 200, width = 400)
    plot(#seq.int(from=5001, to=10000, by = 5),
         mass_parameter,
         type = 'l', xlab = "Iteration", ylab = "Mass parameter")
    dev.off()
    
    pdf(paste0("../figures/mass_parameter_rppa_posterior_chain",
               chain,
               var_sel,
               ".pdf"),
        height = 200, width = 400)
    hist(mass_parameter[1001:2000],
         breaks = 50,
         main = "",
         xlab = "Mass parameter")
    dev.off()
    
    pdf(paste0("../figures/n_clusters_rppa_chain",
               chain,
               var_sel,
               ".pdf"),
        height = 200, width = 400)
    plot(#seq.int(from=5001, to=10000, by = 5),
         n_clusters,
         type = 'l',
         xlab = "Iteration",
         ylab = "Number of clusters")
    dev.off()
    
    pdf(paste0("../figures/n_clusters_rppa_posterior_chain",
               chain,
               var_sel,
               ".pdf"),
        height = 200, width = 400)
    hist(n_clusters[1001:2000],
         # breaks = 1,
         main = "",
         xlab = "Number of clusters")
    dev.off()
  }
}
