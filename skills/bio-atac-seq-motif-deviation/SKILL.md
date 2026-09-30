---
name: bio-atac-seq-motif-deviation
description: Analyze TF motif accessibility variability across samples or single cells using chromVAR. Use when identifying TF motifs whose accessibility correlates with conditions, computing per-sample motif z-scores after matched background correction, comparing to ArchR / Signac equivalents, or distinguishing motif-accessibility signal from per-site footprinting.
license: MIT
category: Data Analysis
author: GPTomics
---

# Motif Deviation (chromVAR)

**"Which TF motifs explain accessibility variation across my samples or cells?"** Compute per-sample (or per-cell) deviation z-scores: how far each TF motif's accessibility departs from expectation, controlling for GC content and overall accessibility via matched background peak sets.

- R: `chromVAR::computeDeviations(counts, motifs)` gives per-sample z-scores.
- R: `chromVAR::computeVariability(dev)` ranks motifs by variability.
- Single-cell: chromVAR called directly on Signac peak counts (`Signac::RunChromVAR()` exists only in Signac 1.16.0 and earlier) or `ArchR::addDeviationsMatrix()`.

chromVAR asks whether peaks containing a motif are systematically more or less accessible than expected. Footprinting (TOBIAS, HINT-ATAC) asks whether one specific motif site is bound. The two are complementary; chromVAR is a summary statistic, so use it when motif site count is well above 100 and use footprinting when individual sites matter.

## What chromVAR computes

For each (motif, sample) pair:
- **Raw deviation**: (observed - expected) / expected, where observed is the count summed over motif-containing peaks and expected is that motif's share of the population-average peak proportions times the sample's depth.
- **Bias-corrected deviation**: raw deviation minus the mean raw deviation of matched-background peak sets (matched on GC content and mean accessibility).
- **Z-score**: bias-corrected deviation divided by the SD of the background raw deviations; the principal output. Positive means the motif is more accessible in this sample than the population average, negative means less. Magnitude depends on the contrast: between GM12878 and K562 the top motifs reached |z| of about 7 per sample.

## Choose a workflow

| Setting | Workflow |
|---------|---------|
| Bulk, 6+ samples, condition contrast | chromVAR, then limma on z-scores; rank by `adj.P.Val` |
| Bulk, 3-5 samples | chromVAR variability ranking only; differential is underpowered |
| Bulk, one condition replicated | chromVAR is uninformative; use footprinting or differential accessibility |
| Bulk time course (5+ points) | z-scores, then spline regression on time |
| scATAC, Signac | `AddMotifs`, chromVAR on the peak counts into a `chromvar` assay, then `FindAllMarkers` (see single-cell reference) |
| scATAC, ArchR | `addPeakMatrix`, `addMotifAnnotations`, `addBgdPeaks`, `addDeviationsMatrix`, `getMarkerFeatures` |
| scATAC + scRNA | chromVAR plus paired DE; consider SCENIC+ for TF-to-target inference |
| Plant / non-model organism | chromVAR with custom PFMs (for example CIS-BP) and a custom BSgenome |

## Bulk workflow

1. Confirm the peakset is at full ATAC scale (5000+ peaks, typically 50k-200k) and record each sample's total library depth (mapped reads or fragments, same unit as the counts) in `depth.tsv`.
2. Build a `SummarizedExperiment` from the peak-by-sample count matrix and peak ranges, set `colData(se)$depth` from `depth.tsv` (total library size, not `colSums(counts)`, which makes the in-peaks fraction 1 and disables that filter); run `addGCBias` with the matching BSgenome (peaks need `seqlengths`).
3. `filterSamples(min_depth=1500, min_in_peaks=0.15)` (total reads per sample >= 1500; fraction of reads in peaks >= 0.15) and `filterPeaks(non_overlapping=TRUE, min_fragments_per_peak=10)`.
4. Fetch motif PFMs (JASPAR vertebrate CORE by default) and run `matchMotifs(..., p.cutoff=5e-5)`. With JASPAR2024, pass the SQLite handle to `getMatrixSet` because TFBSTools does not dispatch on the JASPAR2024 object directly (TFBSTools issue #39).
5. `set.seed()`, then `getBackgroundPeaks(niterations=50)` (background sampling is random: unseeded reruns changed z-scores by up to 1.25 and the significant-motif set by about 12%; report the seed and `niterations`), then `computeDeviations`, `deviationScores` (z-scores; `deviations()` returns raw bias-corrected deviations), and `computeVariability`.
6. For a condition contrast, run limma `lmFit` and `eBayes` on the z-score matrix with a factor design whose first level is the reference, then `topTable(coef=2)`. Use `adj.P.Val` (BH FDR); limma returns no `FDR` column. `logFC` is the z-score difference between groups; `|logFC| >= 0.5` is only a demonstration filter (in GM12878 vs K562, 742 of 879 motifs passed it).

The complete runnable implementation, including variability and heatmap plots and CSV outputs, is [`scripts/chromvar_bulk_analysis.R`](scripts/chromvar_bulk_analysis.R). It reads `peaks.bed`, `counts.tsv` and `depth.tsv` from the working directory and writes motif labels as `ID (name)`; edit the genome and the hard-coded 6-sample `condition` vector before running.

## Guardrails

- Heuristics inherited from the source, not tested or cited here: 5000 peaks, 500 cells per cluster, per-cluster runs above 5x accessibility difference, the 3-5 versus 6+ sample rules, and chromVAR plus scBasset intersection. Only the 6-sample bulk and the two single-cell chains were executed; time-course splines, custom PFMs and scBasset were not.
- Below 5000 peaks or 1500 total reads per sample (bulk), or 500 cells per cluster (single cell), z-scores may be unreliable; all-high variability (above 5) with GC-dominated top motifs signals this. Aggregate cells or enlarge the peakset.
- Keep background defaults (`niterations=50`, `bs=50`) unless benchmarking; do not go below 30 iterations.
- When global accessibility differs more than 5x between cell types, run chromVAR per cluster.
- Run `AddMotifs` and the Signac chromVAR step only after the peakset is final; rerun if peaks change.
- ArchR `getMarkerFeatures` on `MotifMatrix` fails on NA z-scores from very sparse cells; check and drop them first (see single-cell reference).
- Compare z-scores across tools or studies only after recomputing on the same peakset and motif database; z-scores are tool-specific.
- Report the motif database used (JASPAR versus CIS-BP); results differ by database.

## Routed material

- [`references/method-reference.md`](references/method-reference.md): tool taxonomy, chromVAR versus footprinting, per-tool failure modes, variability interpretation, background-matching parameters, tool reconciliation, scBasset alternative, common errors, citations.
- [`references/single-cell.md`](references/single-cell.md): Signac and ArchR code.
- [`references/usage-guide.md`](references/usage-guide.md): prerequisites, inputs, request examples.
- [`scripts/chromvar_bulk_analysis.R`](scripts/chromvar_bulk_analysis.R): bulk chromVAR, variability, limma differential, plots.

## Version compatibility

Executed in the 2026-09-30 audit and fix runs on R 4.4.3 / Bioconductor 3.20: chromVAR 1.30.1, motifmatchr 1.30.0, JASPAR2024 0.99.6, TFBSTools 1.44.0, limma 3.62.1, ArchR 1.0.3, Seurat 5.5.1, Signac 1.17.1 (single-cell reference route) and 1.16.0 (the last release with `RunChromVAR`). Newer Bioconductor releases were not run. With ggplot2 4.x, chromVAR's `plotVariability` warns about deprecated `aes_string()`; the plot still renders. Check `packageVersion('<pkg>')` and `?function_name` before relying on parameters.

## License and provenance

This derived Skill retains GPTomics/bioSkills material by Domen Jemec, from commit `d91ed3d563019e649dc854c56ccd62551359488a` (`atac-seq/motif-deviation`). See [`LICENSE`](LICENSE) for the preserved MIT notice.

## Related Skills

atac-seq/footprinting, atac-seq/differential-accessibility, atac-seq/single-cell-atac, atac-seq/co-accessibility, atac-seq/deep-learning-atac, gene-regulatory-networks/scenic-regulons, chip-seq/motif-analysis, single-cell/clustering.
