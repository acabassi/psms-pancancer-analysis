##################### Pancancer study (Hoadley et al. 2014) ####################
############################ PSMs of copy number data ##########################

rm(list=ls())

library(coda)
library(ComplexHeatmap)
library(klic)
library(mclust)

dataset <- "RPPA"
max_clusters <- 50
n_samples <- 50000
noCuda <- "_noCuda" # Can be either "_noCuda" ("_temp" for mRNA) or ""

for(var_sel in c("", "_alpha1")){
  # Can be a subset of "", "_alpha1", "_alpha01", "_alpha05"
  for(chain in 1:5){

#### Check that maximum number of clusters has been set correctly in mdipp #####

    mcmc_samples <- read.csv(paste0("../mdi/", dataset,"_", max_clusters,
                                    "max_chain", chain, var_sel, "_",
                                    n_samples/1000,
                                    "k", noCuda, ".csv"))
    
####################### Geweke's convergence diagnostic ########################

    mcmc_object <- mcmc(mcmc_samples[(n_samples/10+1):(n_samples/5),1])
    geweke_zscore <- geweke.diag(mcmc_object)$z
    geweke_pval <- 2*pnorm(-abs(geweke_zscore))
    # geweke.plot(mcmc_object)
    
################## Save number of clusters at each iteration ###################

    n_clusters <- rep(0, dim(mcmc_samples)[1])
    for(i in 1:dim(mcmc_samples)[1]){
      sample_i <- mcmc_samples[i,]
      names(sample_i) <- NULL
      sample_i <- unlist(sample_i)
      n_clusters[i] <- length(unique(sample_i))
    }
    plot(n_clusters)

################################## Compute PSM #################################

    mass_parameter <- mcmc_samples[,1]
    mcmc_samples <- mcmc_samples[(n_samples/10+1):(n_samples/5), 2:2422]
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
    save(psm_cn, coph_corr_cn, n_clusters, clusters, ari, geweke_pval,
         file = paste0("../mdi/psm_", dataset, "_chain", chain, var_sel,
                       noCuda, ".RData"))
    # load(paste0("../mdi/psm_", dataset, "_chain", chain, var_sel, noCuda,
    #             ".RData"))

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
    
    pdf(paste0("../figures/psm_mdi_", dataset, "_chain",
               chain, var_sel, noCuda, ".pdf"),
        height = 9,
        width = 8
    )
    draw(
      Hpsm,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom"
    )
    dev.off()

########################### Convergence  assessment ############################

    pdf(paste0("../figures/mass_parameter_", dataset, "_chain",
               chain, var_sel, noCuda, ".pdf"),
        height = 2.5, width = 4)
    plot(mass_parameter[(n_samples/10+1):(n_samples/5)],
         type = 'l', xlab = "Iteration", ylab = "Mass parameter")
    dev.off()
    
    pdf(paste0("../figures/mass_parameter_", dataset,"_posterior_chain",
               chain, var_sel, noCuda, ".pdf"),
        height = 2.5, width = 4)
    hist(mass_parameter[(n_samples/10+1):(n_samples/5)],
         breaks = 50,
         main = "",
         xlab = "Mass parameter")
    dev.off()
    
    pdf(paste0("../figures/n_clusters_", dataset,"_chain",
               chain, var_sel, noCuda, ".pdf"),
        height = 2.5, width = 4)
    plot(n_clusters[(n_samples/10+1):(n_samples/5)],
         type = 'l',
         xlab = "Iteration",
         ylab = "Number of clusters")
    dev.off()
    
    pdf(paste0("../figures/n_clusters_", dataset, "_posterior",
               chain, var_sel, noCuda, ".pdf"),
        height = 2.5, width = 4)
    hist(n_clusters[(n_samples/10+1):(n_samples/5)],
         breaks = 6,
         main = "",
         xlab = "Number of clusters")
    dev.off()
  }
}
