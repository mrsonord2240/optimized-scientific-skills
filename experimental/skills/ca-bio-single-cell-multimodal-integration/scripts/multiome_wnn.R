#!/usr/bin/env Rscript
# Purpose: 10x Multiome (RNA + ATAC, same cell) in Seurat/Signac: RNA PCA, ATAC TF-IDF + LSI, WNN, clusters, UMAP.
# Usage:   Rscript multiome_wnn.R in.rds out_prefix [lsi_start=2] [rna_dims=30] [lsi_end=30]
#            lsi_start=2 drops LSI_1, which usually tracks sequencing depth. Look at out_prefix_depthcor.pdf
#            and pass lsi_start=1 if LSI_1 does not.
# Input:   Seurat object with assays RNA (counts) and ATAC (ChromatinAssay on ONE peak set common to all cells).
# Output:  out_prefix.rds (pca, lsi, weighted.nn, wnn.umap, wsnn clusters), out_prefix_depthcor.pdf,
#          out_prefix_weights.csv (RNA.weight, ATAC.weight per cell).
# Checked: Seurat 5.x, Signac 1.17 (not yet executed after the recut).
suppressPackageStartupMessages({library(Seurat); library(Signac)})

args <- commandArgs(trailingOnly = TRUE)
opts <- list(lsi_start = '2', rna_dims = '30', lsi_end = '30')
positional <- character(0)
for (a in args) {
  kv <- strsplit(a, '=', fixed = TRUE)[[1]]
  if (length(kv) >= 2 && kv[1] %in% names(opts)) opts[[kv[1]]] <- paste(kv[-1], collapse = '=')
  else if (grepl('=', a, fixed = TRUE)) stop('unknown argument ', a)
  else positional <- c(positional, a)
}
if (length(positional) != 2) stop('usage: multiome_wnn.R in.rds out_prefix [lsi_start=2] [rna_dims=30] [lsi_end=30]')
obj <- readRDS(positional[1]); out <- positional[2]
for (a in c('RNA', 'ATAC')) if (!a %in% Assays(obj)) stop('object has no assay ', a)

DefaultAssay(obj) <- 'RNA'
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |> RunPCA()

DefaultAssay(obj) <- 'ATAC'
obj <- RunTFIDF(obj) |> FindTopFeatures(min.cutoff = 'q0') |> RunSVD()
pdf(paste0(out, '_depthcor.pdf'), 6, 4); print(DepthCor(obj)); dev.off()

lsi_dims <- as.integer(opts$lsi_start):as.integer(opts$lsi_end)
obj <- FindMultiModalNeighbors(obj, reduction.list = list('pca', 'lsi'),
                               dims.list = list(1:as.integer(opts$rna_dims), lsi_dims))
obj <- RunUMAP(obj, nn.name = 'weighted.nn', reduction.name = 'wnn.umap')
obj <- FindClusters(obj, graph.name = 'wsnn', algorithm = 3)  # 3 = SLM, the tutorial choice; not Leiden

write.csv(data.frame(cell = colnames(obj), RNA.weight = obj$RNA.weight, ATAC.weight = obj$ATAC.weight),
          paste0(out, '_weights.csv'), row.names = FALSE)
saveRDS(obj, paste0(out, '.rds'))
message('wrote ', out, '.rds (LSI dims ', opts$lsi_start, ':', opts$lsi_end, ')')
