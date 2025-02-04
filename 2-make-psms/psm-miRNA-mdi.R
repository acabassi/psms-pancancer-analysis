##################### Pancancer study (Hoadley et al. 2014) ####################
############################ PSM of copy number data ###########################

rm(list=ls())

library(ComplexHeatmap)
library(klic)
library(mclust)

for(chain in 1:5){
  for(var_sel in c("", "_alpha01", "_alpha05", "_alpha1")){
    
#### Check that maximum number of clusters has been set correctly in mdipp #####

    mcmc_samples <- read.csv(paste0("../mdi/miRNA_100max_chain",
                                    chain,
                                    var_sel,
                                    ".csv"))
    n_clusters <- rep(0, 1000)
    for(i in 1:dim(mcmc_samples)[1]){
      sample_i <-mcmc_samples[i,2:2422]
      names(sample_i) <- NULL
      sample_i <- unlist(sample_i)
      n_clusters[i] <- length(unique(sample_i))
    }
    plot(n_clusters)

################################# Compute PSM ##################################

    mass_parameter <- mcmc_samples[,1]
    mcmc_samples <- mcmc_samples[1001:2000, 2:2422]
    colnames(mcmc_samples) <- sub("Dataset1_", "", colnames(mcmc_samples))
    Rcpp::sourceCpp("makePSM.cpp")
    psm_mirna <- makePSM(as.matrix(mcmc_samples))
    rownames(psm_mirna) <- colnames(psm_mirna) <- colnames(mcmc_samples)
    coph_corr_mirna <- copheneticCorrelation(psm_mirna)

################################## Cluster PSM #################################

    load("../data/samples.RData")
     
    table(n_clusters)
    
    parameters <- list(
      cluster_count =
        as.numeric(names(table(n_clusters))[which.max(table(n_clusters))]))
    psm_mirna <- spectrumShift(psm_mirna, verbose = TRUE, shift = 0.000001)
    kkm_output <- kkmeans(psm_mirna, parameters)
    clusters <- kkm_output$clustering
    names(clusters) <- gsub("\\.", "-", rownames(psm_mirna))
    
    ari <- adjustedRandIndex(clusters, anno_col[names(clusters),]$Tissue)
    save(psm_mirna, coph_corr_mirna, n_clusters, clusters, ari,
         file = paste0("../mdi/psm_mirna_chain", chain, var_sel, ".RData"))
    load(paste0("../mdi/psm_mirna_chain", chain, var_sel, ".RData"))

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
    
    rownames(psm_mirna) <- rownames(psm_mirna) <-
      gsub("\\.", "-", rownames(psm_mirna))
    annotations_mirna <- anno_col[rownames(anno_col) %in% rownames(psm_mirna),]
    annotations_mirna <- annotations_mirna[rownames(psm_mirna),]
    
    Hclusters <-
      rowAnnotation(Tissue = as.character(annotations_mirna$Tissue),
                    name = "Final clusters",
                    annotation_legend_param =
                      my_anno_legend_param("Tissue",nrow=2),
                    show_annotation_name = FALSE,
                    col = list(Tissue =  palette12))
    
    Hpsm <- Heatmap(psm_mirna,
                    col = myBlues,
                    show_column_names = FALSE,
                    show_row_names = FALSE,
                    heatmap_legend_param = my_heatmap_legend_param,
                    name = "Similarity",
                    right_annotation = Hclusters,
                    show_column_dend = FALSE,
                    show_row_dend = FALSE
    )
    
    png(paste0("../figures/psm_mdi_mirna_chain", chain, var_sel, ".png"),
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

    png(paste0("../figures/mass_parameter_mirna_chain",
               chain,
               var_sel,
               ".png"),
        height = 200, width = 400)
    plot(#seq.int(from=5001, to=10000, by = 5),
         mass_parameter[1001:2000],
         type = 'l', xlab = "Iteration", ylab = "Mass parameter")
    dev.off()
    
    png(paste0("../figures/mass_parameter_mirna_posterior_chain",
               chain,
               var_sel,
               ".png"),
        height = 200, width = 400)
    hist(mass_parameter[1001:2000],
         breaks = 50,
         main = "",
         xlab = "Mass parameter")
    dev.off()
    
    png(paste0("../figures/n_clusters_mirna_chain",
               chain,
               var_sel,
               ".png"),
        height = 200, width = 400)
    plot(#seq.int(from=5001, to=10000, by = 5),
         n_clusters,
         type = 'l',
         xlab = "Iteration",
         ylab = "Number of clusters")
    dev.off()
    
    png(paste0("../figures/n_clusters_mirna_posterior_chain",
               chain,
               var_sel,
               ".png"),
        height = 200, width = 400)
    hist(n_clusters[1001:2000],
         breaks = 6,
         main = "",
         xlab = "Number of clusters")
    dev.off()
  }
}
