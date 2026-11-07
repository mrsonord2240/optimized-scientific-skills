#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop('usage: cluster_peaks.R clustered_signac.rds /path/to/macs3 output.rds')
if (!file.exists(args[1]) || !file.exists(args[2]) || file.exists(args[3]))
    stop('object and MACS3 must exist; output must be fresh')
suppressPackageStartupMessages({ library(Signac); library(Seurat) })
obj <- readRDS(args[1])
if (!'seurat_clusters' %in% colnames(obj[[]])) stop('object needs seurat_clusters metadata')
peaks <- CallPeaks(obj, group.by='seurat_clusters', macs2.path=normalizePath(args[2]),
                   cleanup=FALSE, format='BED', shift=-75, extsize=150,
                   additional.args='-p 0.01')
saveRDS(peaks, args[3])
