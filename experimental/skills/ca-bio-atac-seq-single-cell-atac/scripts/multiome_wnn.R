#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop('usage: multiome_wnn.R paired_multiome.rds output.rds')
if (!file.exists(args[1]) || file.exists(args[2])) stop('input must exist and output must be fresh')
suppressPackageStartupMessages({ library(Signac); library(Seurat) })
obj <- readRDS(args[1])
if (!all(c('RNA', 'ATAC') %in% names(obj))) stop('object must contain paired RNA and ATAC assays')
DefaultAssay(obj) <- 'RNA'
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
DefaultAssay(obj) <- 'ATAC'
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff='q0')
obj <- RunSVD(obj)
obj <- FindMultiModalNeighbors(obj, reduction.list=list('pca', 'lsi'), dims.list=list(1:30, 2:30))
obj <- RunUMAP(obj, nn.name='weighted.nn', reduction.name='wnn.umap')
saveRDS(obj, args[2])
