# Seurat v5 bridge integration: map an unpaired scATAC query onto a labelled scRNA reference
# through a paired multiome bridge dataset.
# Inputs (Seurat objects saved with saveRDS):
#   rna.rds   labelled scRNA reference: meta.data column 'celltype'; NormalizeData/ScaleData/
#             RunPCA ('pca'); RunUMAP(return.model = TRUE) ('umap')
#   multi.rds paired multiome bridge with 'RNA' (normalized) and 'ATAC' (RunTFIDF + RunSVD -> 'lsi')
#   atac.rds  unpaired scATAC query on the SAME peak set as the bridge's ATAC assay, RunTFIDF only
# Output: <out.rds> = the query with predicted.celltype, predicted.celltype.score and ref.umap
# Usage:  Rscript scripts/seurat_bridge_integration.R rna.rds multi.rds atac.rds out.rds [ndims] [first_lsi_dim] [SCT|LogNormalize]
# If Seurat/Signac live in a private library, set BIO_SKILLS_R_LIB to that directory before invoking Rscript.
# Checked on Seurat 5.5.0, Signac 1.17.1.
private_lib <- Sys.getenv('BIO_SKILLS_R_LIB', unset = '')
if (nzchar(private_lib)) {
  .libPaths(c(strsplit(private_lib, .Platform$path.sep, fixed = TRUE)[[1]], .libPaths()))
}
required_packages <- c('Seurat', 'Signac')
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop('Missing required R package(s): ', paste(missing_packages, collapse = ', '),
       '. Install them in an R library, or set BIO_SKILLS_R_LIB to that library before running this script.',
       call. = FALSE)
}
suppressMessages({library(Seurat); library(Signac)})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 4)
ndims <- if (length(args) >= 5) as.integer(args[5]) else 20L
lsi_start <- if (length(args) >= 6) as.integer(args[6]) else 2L     # 2 drops LSI_1: ONLY if DepthCor confirms it tracks depth
norm_method <- if (length(args) >= 7) args[7] else 'LogNormalize'   # 'SCT' if the reference and bridge RNA were SCTransformed
rna <- readRDS(args[1]); multi <- readRDS(args[2]); atac <- readRDS(args[3])

bridge <- PrepareBridgeReference(
    reference = rna, bridge = multi,
    reference.reduction = 'pca', reference.dims = 1:ndims,
    normalization.method = norm_method,
    bridge.ref.assay = 'RNA', bridge.query.assay = 'ATAC',
    supervised.reduction = 'slsi', laplacian.reduction.dims = 1:ndims)

anchors <- FindBridgeTransferAnchors(extended.reference = bridge, query = atac,
                                     reduction = 'lsiproject', dims = lsi_start:ndims)

# reference = the object PrepareBridgeReference RETURNED, not the original scRNA object
# (the anchorset lives in its 'Bridge' assay; passing the original errors "assay ... does not match")
atac <- MapQuery(anchorset = anchors, reference = bridge, query = atac,
                 refdata = list(celltype = 'celltype'), reduction.model = 'umap')

saveRDS(atac, args[4])
cat('Mapped', ncol(atac), 'query cells; labels:', paste(names(table(atac$predicted.celltype)), collapse = ', '), '\n')
