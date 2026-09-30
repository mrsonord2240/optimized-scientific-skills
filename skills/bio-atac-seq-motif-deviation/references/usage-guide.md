# Motif Deviation (chromVAR) usage guide

## Prerequisites

```r
BiocManager::install(c('chromVAR', 'motifmatchr', 'JASPAR2024', 'TFBSTools',
                       'BSgenome.Hsapiens.UCSC.hg38', 'SummarizedExperiment',
                       'limma', 'Signac', 'Seurat'))

# ArchR alternative (if using ArchR ecosystem)
remotes::install_github('GreenleafLab/ArchR', ref='master', repos=BiocManager::repositories())
```

Inputs: a peak count matrix (rows = peaks, columns = samples) or a Seurat/ArchR object; peak ranges (BED or GRanges); a motif PFM database (JASPAR 2024 default; HOCOMOCO or CIS-BP for special cases). The bulk script also uses pheatmap and RSQLite.

## Example requests

- "Run chromVAR on this peak count matrix with JASPAR 2024 vertebrate CORE motifs. Add GC bias correction, filter samples below 1500 total reads or with < 0.15 of reads in peaks (I will supply `depth.tsv` with total reads per sample), then report the top 20 most variable motifs and their z-score matrix."
- "Run limma on the chromVAR z-score matrix to identify motifs differing between control and treated. Report adj.P.Val < 0.05 (and abs(logFC) > 0.5 as a demonstration filter); do not use raw counts."
- "Add motifs to my Seurat object after the peakset is finalized, run chromVAR on the peak counts (`RunChromVAR` only exists in Signac <= 1.16), then FindAllMarkers on the chromvar assay with `mean.fxn=rowMeans` and `fc.name='avg_diff'`."
- "In my ArchR project, add reproducible peakset, addPeakMatrix, addMotifAnnotations with `motifSet='cisbp'`, addBgdPeaks, addDeviationsMatrix; then getMarkerFeatures on MotifMatrix grouped by Clusters."
- "My chromVAR variability scores are all above 5; verify peak count is at full ATAC scale (50k+) and per-sample total reads exceed 1500. If too sparse, aggregate cells before running."
- "Fit a spline regression on chromVAR z-scores across the time course; identify motifs with non-monotonic trajectories."

## Expected steps

Verify peakset and depth; build the SummarizedExperiment; add GC bias; filter samples and peaks; match motifs; sample background peaks; compute deviations and variability; optionally test differential motifs (limma) or per-cluster markers (FindAllMarkers / getMarkerFeatures); plot variability, top-motif heatmap, and PCA.

## Tips

- Variability needs no condition labels and suits unsupervised discovery; differential testing needs labels.
- ArchR `cisbp` gave 870 motifs in the ArchR 1.0.3 test; JASPAR2024 CORE vertebrates gives 879 (1,912 with `all_versions=TRUE`). Choose one database and report it.
- Use limma `adj.P.Val` for differential motifs; raw t-tests on z-scores are mis-calibrated.
- For differential peak counts use DiffBind or differential accessibility, not chromVAR.
