##################### Pancancer study (Hoadley et al. 2014) ####################
############################ PSM of copy number data ###########################

rm(list=ls())

library(ComplexHeatmap)
library(klic)
library(mclust)

for(chain in 1:5){
  for(var_sel in c("_alpha01", "_alpha05", "_alpha1")){ # Also "", 
    
#### Check that maximum number of clusters has been set correctly in mdipp #####

    mcmc_samples <- read.csv(paste0("../mdi/mRNA_250max_chain",
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
    psm_mrna <- makePSM(as.matrix(mcmc_samples))
    rownames(psm_mrna) <- colnames(psm_mrna) <- colnames(mcmc_samples)
    coph_corr_mrna <- copheneticCorrelation(psm_mrna)

################################## Cluster PSM #################################

    load("../data/samples.RData")
    
    table(n_clusters)
    
    parameters <- list(
      cluster_count =
        as.numeric(names(table(n_clusters))[which.max(table(n_clusters))]))
    psm_mrna <- spectrumShift(psm_mrna, verbose = TRUE, shift = 0.000001)
    kkm_output <- kkmeans(psm_mrna, parameters)
    clusters <- kkm_output$clustering
    names(clusters) <- gsub("\\.", "-", rownames(psm_mrna))
    
    ari <- adjustedRandIndex(clusters, anno_col[names(clusters),]$Tissue)
    save(psm_mrna, coph_corr_mrna, n_clusters, clusters, ari,
         file = paste0("../mdi/psm_mRNA_chain", chain, ".RData"))
    load(paste0("../mdi/psm_mRNA_chain", chain, ".RData"))

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
                    annotation_legend_param =
                      my_anno_legend_param("Tissue",nrow=2),
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
    
    pdf(paste0("../figures/psm_mdi_mrna_chain",chain,".pdf"),
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

    pdf(paste0("../figures/mass_parameter_mrna_chain",
               chain,
               var_sel,
               ".pdf"),
        height = 200, width = 400)
    plot(seq.int(from=5001, to=10000, by = 5),
         mass_parameter[1001:2000],
         type = 'l', xlab = "Iteration", ylab = "Mass parameter")
    dev.off()
    
    pdf(paste0("../figures/mass_parameter_mrna_posterior_chain",
               chain,
               var_sel,
               ".pdf"),
        height = 200, width = 400)
    hist(mass_parameter[1001:2000],
         breaks = 50,
         main = "",
         xlab = "Mass parameter")
    dev.off()
    
    pdf(paste0("../figures/n_clusters_mrna_chain",chain,".pdf"),
        height = 200, width = 400)
    plot(seq.int(from=5001, to=10000, by = 5),
         n_clusters,
         type = 'l',
         xlab = "Iteration",
         ylab = "Number of clusters")
    dev.off()
    
    pdf("../figures/n_clusters_mrna_posterior.pdf",
        height = 200, width = 400)
    hist(n_clusters,
         # breaks = 50,
         main = "",
         xlab = "Number of clusters")
    dev.off()
  }
}