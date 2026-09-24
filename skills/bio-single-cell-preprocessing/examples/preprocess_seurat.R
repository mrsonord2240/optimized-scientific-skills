# Reference: Seurat 5.0+, scran 1.30+ | Verify API if version differs
# Preprocess single-cell data with Seurat

library(Seurat)

counts <- Read10X(data.dir = 'filtered_feature_bc_matrix/')
seurat_obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
cat('Raw:', ncol(seurat_obj), 'cells,', nrow(seurat_obj), 'genes\n')

seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = '^MT-')

is_outlier <- function(x, nmads) {
    centre <- median(x)
    spread <- mad(x, constant = 1)
    cat(sprintf('median=%.3f, MAD=%.3f\n', centre, spread))
    if (spread == 0) {
        stop('MAD collapsed; use documented fixed cutoffs instead of MAD filtering')
    }
    x < centre - nmads * spread | x > centre + nmads * spread
}

mito_hard_caps <- c(nuclei = NA_real_, pbmc = 8, cardiac = 30, hepatic = 30,
                    skeletal_muscle = 40, unknown = NA_real_)
tissue <- 'unknown'  # set from sample metadata before filtering
mito_hard_cap <- unname(mito_hard_caps[tissue])
hard_mito <- if (is.na(mito_hard_cap)) rep(FALSE, ncol(seurat_obj)) else seurat_obj$percent.mt > mito_hard_cap

seurat_obj$qc_outlier <- (
    is_outlier(log1p(seurat_obj$nCount_RNA), 5) |
    is_outlier(log1p(seurat_obj$nFeature_RNA), 5) |
    is_outlier(seurat_obj$percent.mt, 3) |
    hard_mito
)
survival_fraction <- mean(!seurat_obj$qc_outlier)
if (survival_fraction < 0.80) {
    stop(sprintf('Only %.1f%% of barcodes survive QC; inspect MADs and use fixed, tissue-aware cutoffs',
                 100 * survival_fraction))
}
seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[!seurat_obj$qc_outlier])
cat('Filtered:', ncol(seurat_obj), 'cells\n')

# Do not reflexively regress out percent.mt: it is confounded with real cell state and erases biology
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
cat('HVGs:', length(VariableFeatures(seurat_obj)), '\n')

saveRDS(seurat_obj, file = 'preprocessed.rds')
cat('Saved preprocessed data\n')
