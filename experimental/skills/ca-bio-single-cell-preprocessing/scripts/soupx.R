#!/usr/bin/env Rscript
# Purpose: SoupX ambient-RNA removal from a Cell Ranger outs/ directory (raw AND filtered matrices).
# Usage:   Rscript soupx.R cellranger_outs/ adjusted_counts.rds
# Output:  adjusted_counts.rds (integer-rounded dgCMatrix, genes x cells); prints the contamination fraction rho.
# Notes:   load10X imports clusters only when Cell Ranger's analysis output is present; otherwise a quick Seurat
#          pass on the filtered matrix supplies them. autoEstCont refuses to run unclustered.
# Checked: SoupX 1.6.2, Seurat 5.5.0 (R 4.4.3); not yet executed after the recut.
suppressPackageStartupMessages({ library(SoupX) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('usage: soupx.R cellranger_outs/ adjusted_counts.rds')
sc <- load10X(args[1])  # needs BOTH raw and filtered
if (is.null(sc$metaData$clusters)) {
  so <- Seurat::CreateSeuratObject(sc$toc)
  so <- Seurat::NormalizeData(so, verbose = FALSE)
  so <- Seurat::FindVariableFeatures(so, verbose = FALSE)
  so <- Seurat::ScaleData(so, verbose = FALSE)
  so <- Seurat::RunPCA(so, npcs = 30, verbose = FALSE)
  so <- Seurat::FindNeighbors(so, dims = 1:20, verbose = FALSE)
  so <- Seurat::FindClusters(so, resolution = 0.8, verbose = FALSE)
  sc <- setClusters(sc, setNames(as.character(Seurat::Idents(so)), colnames(so)))
}
sc <- autoEstCont(sc)
cat('contamination fraction rho:', sc$fit$rhoEst, '\n')
saveRDS(adjustCounts(sc, roundToInt = TRUE), args[2])
cat('wrote', args[2], '\n')
