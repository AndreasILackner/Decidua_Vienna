### Last updated 1/3/2025

#if (!requireNamespace("BiocManager", quietly=TRUE)) {
#  install.packages("BiocManager")
#}
#BiocManager::install("zellkonverter")
#install.packages('rliger')
#install.packages('RcppPlanc', repos = c('https://welch-lab.r-universe.dev', 'https://cloud.r-project.org'))
#install.packages('Matrix')

library('zellkonverter')
library('rliger')
library('Matrix')

sce1 <- readH5AD("./adata_cite.h5ad", layers = TRUE)
sce2 <- readH5AD("./adata_xen.h5ad", layers = TRUE)

lig <- createLiger(list(citeseq = sce1, xenium = sce2),
                   modal = c("default", "default"))

lig <- normalize(lig)
lig <- selectGenes(lig, thresh = 0.3, combine = 'union', useUnsharedDatasets = 'citeseq')
lig <- scaleNotCenter(lig)
lig <- runUINMF(lig, k = 50, lambda = 10)
lig <- quantileNorm(lig, reference = "citeseq")
lig <- runUMAP(lig)
plotDatasetDimRed(lig)
lig <- runCluster(lig, nNeighbors = 30, saveSNN = TRUE)

norm_factors <- lig@H.norm
write.csv(norm_factors, './norm_factors.csv')

SNN <- lig@uns$snn
writeMM(SNN, file = "./SNN.mtx")

umap_coords <- lig@dimReds$UMAP
write.csv(umap_coords, './umap_coords.csv')
