#!/usr/bin/env Rscript
# Purpose: load 10x Cell Ranger output into a Seurat object and save it as RDS.
# Usage:   Rscript load_10x.R <matrix_dir | file.h5> out.rds [project=PBMC]
# Output:  out.rds. With several feature types (Gene Expression + Antibody Capture ...) the Gene Expression block is
#          the RNA assay and each other block becomes an extra assay named after its feature type.
# Note:    min.cells=3 and min.features=200 drop sparse genes and near-empty barcodes; load the raw matrix with
#          another tool if the empty droplets are needed.
# Checked: Seurat 5.5.0 (R 4.4.3); not yet executed after the recut.
suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('usage: load_10x.R <matrix_dir | file.h5> out.rds [project=PBMC]')
project <- 'PBMC'
for (a in args[-(1:2)]) {
  kv <- strsplit(a, '=', fixed = TRUE)[[1]]
  if (length(kv) != 2 || kv[1] != 'project') stop('unknown argument ', a)
  project <- kv[2]
}
counts <- if (grepl('\\.h5$', args[1])) Read10X_h5(args[1]) else Read10X(data.dir = args[1])  # list when several feature types
if (is.list(counts)) {
  obj <- CreateSeuratObject(counts = counts[['Gene Expression']], project = project, min.cells = 3, min.features = 200)
  for (name in setdiff(names(counts), 'Gene Expression')) {
    extra <- counts[[name]][, colnames(obj), drop = FALSE]
    obj[[gsub(' ', '', name)]] <- CreateAssay5Object(counts = extra)
  }
} else {
  obj <- CreateSeuratObject(counts = counts, project = project, min.cells = 3, min.features = 200)
}
print(obj)
cat('Cells:', ncol(obj), ' Genes:', nrow(obj), '\n')
saveRDS(obj, file = args[2])
cat('wrote', args[2], '\n')
