#!/usr/bin/env Rscript
# Purpose: preprocess single-cell counts with Seurat: MAD QC, normalization, variable features.
# Usage:   Rscript preprocess_seurat.R <10x_dir_or_.h5> out.rds [tissue=unknown] [mt_pattern=^MT-] [normalization=sct]
#            tissue: nuclei|unknown = MAD only; pbmc 8; cardiac 30; hepatic 30; skeletal_muscle 40 (percent.mt cap)
#            normalization: sct = SCTransform; lognorm = NormalizeData(scale.factor 10000) + FindVariableFeatures
# Output:  out.rds. Stops when percent.mt is all zero (symbol pattern against Ensembl IDs), a MAD is 0, or fewer
#          than 80% of cells survive.
# Checked: Seurat 5.5.0 (R 4.4.3); not yet executed after the recut.
suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('usage: preprocess_seurat.R <10x_dir_or_.h5> out.rds [tissue=] [mt_pattern=] [normalization=]')
opts <- list(tissue = 'unknown', mt_pattern = '^MT-', normalization = 'sct')
for (a in args[-(1:2)]) {
  kv <- strsplit(a, '=', fixed = TRUE)[[1]]
  if (length(kv) != 2 || !kv[1] %in% names(opts)) stop('unknown argument ', a)
  opts[[kv[1]]] <- kv[2]
}
caps <- c(nuclei = NA_real_, pbmc = 8, cardiac = 30, hepatic = 30, skeletal_muscle = 40, unknown = NA_real_)
if (!opts$tissue %in% names(caps)) stop('tissue must be one of: ', paste(names(caps), collapse = ', '))
if (!opts$normalization %in% c('sct', 'lognorm')) stop('normalization must be sct or lognorm')
mito_cap <- unname(caps[opts$tissue])

counts <- if (grepl('\\.h5$', args[1])) Read10X_h5(args[1]) else Read10X(data.dir = args[1])
if (is.list(counts)) counts <- counts[['Gene Expression']]
obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
cat('Raw:', ncol(obj), 'cells,', nrow(obj), 'genes\n')

obj[['percent.mt']] <- PercentageFeatureSet(obj, pattern = opts$mt_pattern)
if (all(obj$percent.mt == 0)) stop('percent.mt is 0 for every cell: pattern ', opts$mt_pattern,
                                   ' matches no gene (Ensembl IDs? mouse needs ^mt-)')
is_outlier <- function(x, nmads) {
  centre <- median(x)
  spread <- mad(x, constant = 1)
  cat(sprintf('median=%.3f, MAD=%.3f\n', centre, spread))
  if (spread == 0) stop('MAD collapsed; use documented fixed cutoffs instead of MAD filtering')
  x < centre - nmads * spread | x > centre + nmads * spread
}
hard_mito <- if (is.na(mito_cap)) rep(FALSE, ncol(obj)) else obj$percent.mt > mito_cap
obj$qc_outlier <- (is_outlier(log1p(obj$nCount_RNA), 5) | is_outlier(log1p(obj$nFeature_RNA), 5) |
                   is_outlier(obj$percent.mt, 3) | hard_mito)
survival <- mean(!obj$qc_outlier)
if (survival < 0.80) stop(sprintf('Only %.1f%% of barcodes survive QC; use fixed, tissue-aware cutoffs', 100 * survival))
obj <- subset(obj, cells = colnames(obj)[!obj$qc_outlier])
cat('Filtered:', ncol(obj), 'cells\n')

# Do not regress out percent.mt: it is confounded with real cell state.
if (opts$normalization == 'sct') {
  obj <- SCTransform(obj, verbose = FALSE)
} else {
  obj <- NormalizeData(obj, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
}
cat('HVGs:', length(VariableFeatures(obj)), '\n')
saveRDS(obj, file = args[2])
cat('wrote', args[2], '\n')
