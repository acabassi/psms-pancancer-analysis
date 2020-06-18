############### Select tumour samples to be used in the analysis ###############

rm(list=ls())

load("data/samples_CN.RData")
load("data/samples_mRNA.RData")
load("data/samples_miRNA.RData")
load("data/samples_methylation.RData")
load("data/samples_RPPA.RData")

samples_intersection <-
  intersect(
    intersect(
      intersect(
        intersect(samples_CN, samples_methylation),
        samples_miRNA),
      samples_mRNA),
    samples_RPPA)

save(samples_intersection, file = "data/samples_intersection.RData")
