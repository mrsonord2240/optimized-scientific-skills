#!/usr/bin/env Rscript
# Purpose: cluster a preprocessed Seurat object: PCA, SNN graph, a resolution sweep, UMAP, clustree plot.
# Usage:   Rscript cluster_seurat.R in.rds out.rds [dims=30] [resolutions=0.2,0.4,0.6,0.8,1.0] [algorithm=4] [seed=0]
#            algorithm=4 Leiden (needs leidenbase), 1 Louvain, 3 SLM
# Output:  out.rds (columns RNA_snn_res.<r>), clustree.pdf, <out>.sweep.tsv (cluster count per resolution).
# Input:   normalized object with variable features set (NormalizeData or SCTransform already run).
# Checked: Seurat 5.5.0 (R 4.4.3); not yet executed after the recut.
suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('usage: cluster_seurat.R in.rds out.rds [dims=30] [resolutions=...] [algorithm=4] [seed=0]')
opts <- list(dims = '30', resolutions = '0.2,0.4,0.6,0.8,1.0', algorithm = '4', seed = '0')
for (a in args[-(1:2)]) {
  kv <- strsplit(a, '=', fixed = TRUE)[[1]]
  if (length(kv) != 2 || !kv[1] %in% names(opts)) stop('unknown argument ', a)
  opts[[kv[1]]] <- kv[2]
}
dims <- seq_len(as.integer(opts$dims))
resolutions <- as.numeric(strsplit(opts$resolutions, ',')[[1]])
algorithm <- as.integer(opts$algorithm)
if (algorithm == 4 && !requireNamespace('leidenbase', quietly = TRUE))
  stop("algorithm=4 (Leiden) needs install.packages('leidenbase'); or pass algorithm=1 for Louvain")

obj <- readRDS(args[1])
if (length(VariableFeatures(obj)) == 0) stop('no variable features set; run the preprocessing Skill first')
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
obj <- FindNeighbors(obj, dims = dims, verbose = FALSE)
obj <- FindClusters(obj, resolution = resolutions, algorithm = algorithm, random.seed = as.integer(opts$seed), verbose = FALSE)
obj <- RunUMAP(obj, dims = dims, seed.use = as.integer(opts$seed), verbose = FALSE)

prefix <- paste0(DefaultAssay(obj), '_snn_res.')
cols <- paste0(prefix, resolutions)
cols <- cols[cols %in% colnames(obj@meta.data)]
sweep <- data.frame(resolution = resolutions[paste0(prefix, resolutions) %in% cols], column = cols,
                    n_clusters = sapply(cols, function(c) length(unique(obj@meta.data[[c]]))))
print(sweep)
write.table(sweep, sub('\\.[^.]*$', '.sweep.tsv', args[2]), sep = '\t', quote = FALSE, row.names = FALSE)
if (requireNamespace('clustree', quietly = TRUE)) {
  pdf('clustree.pdf'); print(clustree::clustree(obj, prefix = prefix)); dev.off()
} else message('clustree not installed; skipping clustree.pdf')
saveRDS(obj, args[2])
cat('wrote', args[2], '\n')
