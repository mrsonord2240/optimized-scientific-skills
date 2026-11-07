#!/usr/bin/env Rscript
# Purpose: per-cell chromVAR z-scores for a Signac/Seurat scATAC object, stored as a `chromvar` assay, plus
#   per-cluster differential motifs. Does not use Signac::RunChromVAR (removed in Signac 1.17).
# Usage:   Rscript chromvar_signac.R in.rds out.rds markers.csv [assay=ATAC] [genome=BSgenome.Hsapiens.UCSC.hg38] [seed=2024]
# Input:   Seurat object with a ChromatinAssay on the FINAL peakset and an identity (clusters) set. AddMotifs is run here
#          with JASPAR2024 vertebrates CORE; rerun this script if the peaks change.
# Output:  out.rds with assay `chromvar` (z-scores, motifs x cells); markers.csv = FindAllMarkers on that assay
#          (mean.fxn=rowMeans, fc.name=avg_diff) with the motif name added.
# Checked: Signac 1.17.1, Seurat 5.5.1, chromVAR 1.30.1 (not yet executed after the recut).
suppressPackageStartupMessages({
  library(Signac); library(Seurat); library(JASPAR2024); library(TFBSTools); library(RSQLite)
  library(chromVAR); library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
opts <- list(assay = 'ATAC', genome = 'BSgenome.Hsapiens.UCSC.hg38', seed = '2024')
positional <- character(0)
for (a in args) {
  kv <- strsplit(a, '=', fixed = TRUE)[[1]]
  if (length(kv) >= 2 && kv[1] %in% names(opts)) opts[[kv[1]]] <- paste(kv[-1], collapse = '=')
  else if (grepl('=', a, fixed = TRUE)) stop('unknown argument ', a)
  else positional <- c(positional, a)
}
if (length(positional) != 3) stop('usage: chromvar_signac.R in.rds out.rds markers.csv [assay=ATAC] [genome=...] [seed=2024]')
if (!requireNamespace(opts$genome, quietly = TRUE)) stop('Install ', opts$genome, ' first')
genome <- getExportedValue(opts$genome, opts$genome)

seurat_obj <- readRDS(positional[1])
jaspar2024 <- JASPAR2024::JASPAR2024()
sq <- dbConnect(SQLite(), db(jaspar2024))  # TFBSTools issue #39
pfm <- getMatrixSet(sq, opts = list(collection = 'CORE', tax_group = 'vertebrates'))
seurat_obj <- AddMotifs(seurat_obj, genome = genome, pfm = pfm, assay = opts$assay)

counts <- GetAssayData(seurat_obj, assay = opts$assay, layer = 'counts')
se <- SummarizedExperiment(assays = list(counts = counts), rowRanges = granges(seurat_obj[[opts$assay]]))
se <- addGCBias(se, genome = genome)
keep <- Matrix::rowSums(counts) > 0 & !is.na(rowData(se)$bias)  # empty peaks and undefined GC break background matching
se <- se[keep, ]
set.seed(as.integer(opts$seed))
bg <- getBackgroundPeaks(se)
dev <- computeDeviations(se, annotations = GetMotifData(seurat_obj[[opts$assay]], slot = 'data')[keep, ], background_peaks = bg)
z <- deviationScores(dev)
cat('NA z-scores:', sum(is.na(z)), '\n')
if (anyNA(z)) warning('NA z-scores present: sparse cells can give 0/0; FindAllMarkers may fail. Count and drop those cells and report how many.')
seurat_obj[['chromvar']] <- CreateAssayObject(data = z)
DefaultAssay(seurat_obj) <- 'chromvar'

markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, mean.fxn = rowMeans, fc.name = 'avg_diff')
markers$motif_name <- ConvertMotifID(seurat_obj[[opts$assay]], id = markers$gene)
write.csv(markers, positional[3], row.names = FALSE)
saveRDS(seurat_obj, positional[2])
message('wrote ', positional[2], ' and ', positional[3], ' (seed ', opts$seed, ')')
