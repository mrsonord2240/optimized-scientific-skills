# Single-cell chromVAR code

Tested on Signac 1.17.1 and ArchR 1.0.3 (10x PBMC 5k, chr1:1-30 Mb slice).

## Signac

**Goal:** Compute per-cell TF-motif z-scores in a Seurat scATAC workflow and call cluster-marker motifs.

**Approach:** Attach motifs with `AddMotifs`, run chromVAR directly on the peak counts with the Signac motif matrix, store the z-scores as a `chromvar` assay, then `FindAllMarkers` with `mean.fxn=rowMeans`. `Signac::RunChromVAR()` was removed in Signac 1.17; it remains available in 1.16.0 and earlier.

```r
library(Signac); library(Seurat); library(JASPAR2024); library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38); library(RSQLite)
library(chromVAR); library(SummarizedExperiment)

# Assume `seurat_obj` has an ATAC assay with consensus peaks
# JASPAR2024 + TFBSTools workaround (see TFBSTools issue #39):
jaspar2024 <- JASPAR2024::JASPAR2024()
sq <- dbConnect(SQLite(), db(jaspar2024))
pfm <- getMatrixSet(sq, opts=list(collection='CORE', tax_group='vertebrates'))
seurat_obj <- AddMotifs(seurat_obj, genome=BSgenome.Hsapiens.UCSC.hg38, pfm=pfm)

# chromVAR on the ATAC counts
counts <- GetAssayData(seurat_obj, assay='ATAC', layer='counts')
se <- SummarizedExperiment(assays=list(counts=counts), rowRanges=granges(seurat_obj[['ATAC']]))
se <- addGCBias(se, genome=BSgenome.Hsapiens.UCSC.hg38)
# empty peaks and peaks with undefined GC (assembly gaps, NaN bias) break background matching
keep <- Matrix::rowSums(counts) > 0 & !is.na(rowData(se)$bias)
se <- se[keep, ]
set.seed(2024)                                   # background sampling is random; report the seed
bg <- getBackgroundPeaks(se)
dev <- computeDeviations(se, annotations=GetMotifData(seurat_obj[['ATAC']], slot='data')[keep, ],
                         background_peaks=bg)
z <- deviationScores(dev)                        # motifs x cells
cat('NA z-scores:', sum(is.na(z)), '\n')       # expect 0; see the ArchR note if not
seurat_obj[['chromvar']] <- CreateAssayObject(data=z)
DefaultAssay(seurat_obj) <- 'chromvar'

# Per-cluster differential motifs.
# `mean.fxn` is the standard FindAllMarkers/FindMarkers control for the per-feature summary.
# `fc.name` controls the output column name and is accepted by Seurat 4.x/5.x; if it errors,
# fall back to renaming the output column post-hoc.
markers <- FindAllMarkers(seurat_obj, only.pos=TRUE, mean.fxn=rowMeans, fc.name='avg_diff')
```

`mean.fxn=rowMeans` is required for z-score-style data; the default fold-change function (designed for log-counts) is not meaningful on chromVAR z-scores. Rows of the `chromvar` assay are JASPAR IDs; map them with `ConvertMotifID(seurat_obj[['ATAC']], id=markers$gene)`.

## ArchR

**Goal:** Compute per-cell TF-motif deviations and per-cluster marker motifs within ArchR.

**Approach:** Build the reproducible peakset, attach CIS-BP motif annotations, sample matched background peaks, run `addDeviationsMatrix`, and call `getMarkerFeatures` on the MotifMatrix per cluster.

```r
library(ArchR)
proj <- addReproduciblePeakSet(proj, groupBy='Clusters', pathToMacs2='/path/macs2')
proj <- addPeakMatrix(proj)
proj <- addMotifAnnotations(proj, motifSet='cisbp', name='Motif')
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(proj, peakAnnotation='Motif')

# Sparse cells x rare motifs can give 0/0 (NaN) z-scores, and getMarkerFeatures (presto) refuses NA.
mm <- getMatrixFromProject(proj, useMatrix='MotifMatrix')
na_cells <- colnames(mm)[colSums(is.na(assays(mm)$z)) > 0]
length(na_cells)
if (length(na_cells) > 0) proj <- subsetCells(proj, cellNames=setdiff(proj$cellNames, na_cells))

# Per-cluster deviation summary
markersMotifs <- getMarkerFeatures(proj, useMatrix='MotifMatrix',
                                   groupBy='Clusters', useSeqnames='z')
```

ArchR uses `cisbp` by default (870 motifs in the ArchR 1.0.3 test); `JASPAR2020` is the alternative set. Use `matrixName='PeakMatrix'` after `addReproduciblePeakSet`; tile-based deviations add intergenic noise and are mainly for embedding.

The NA cells in the slice test had 26-203 reads in the peak matrix (median cell 917); their rare motifs had 3-12 matched peaks with no reads, so observed and background deviations were identical and the background SD was zero. Dropping them let `getMarkerFeatures` finish. Whether this recurs on a full-depth dataset is untested; count NA first and report how many cells were dropped.
