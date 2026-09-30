# Motif deviation method reference

Methodology evolves; verify against current chromVAR (Schep 2017), ArchR (Granja 2021), and Signac (Stuart 2021) benchmarks before locking pipelines.

## Algorithmic taxonomy

| Tool | Input | Background | Output | Best for | Fails when |
|------|-------|------------|--------|----------|------------|
| chromVAR | Peak count matrix + motif annotations | Matched GC + accessibility (50 peaks per match by default) | Per-sample motif z-score | Bulk + single-cell (sparse-aware); cross-sample variability | < 1500 total reads/sample (bulk) or < 500 cells/cluster (sc); too few peaks (< 5000) |
| Signac + chromVAR (`RunChromVAR` in Signac <= 1.16) | Seurat scATAC object | Same as chromVAR | Motif assay in Seurat object | Seurat-ecosystem workflows | Same as chromVAR; needs Seurat object setup; `RunChromVAR` absent in Signac >= 1.17 |
| ArchR::addDeviationsMatrix | ArrowFile + tile/peak matrix | ArchR getBgdPeaks (matched on GC + log accessibility) | Per-cell deviation matrix in ArchR project | ArchR ecosystem; faster on large scATAC | ArchR-specific format; not portable to chromVAR objects |
| Signac::FindMarkers (motifs as features) | Motif accessibility matrix (chromvar assay) | Per-cell-cluster | Differential motifs per cluster | Cluster-level differential | Test must be on z-scores; raw counts mislead |
| SCENIC+ | Expression + accessibility | Multi-modal | TF-to-target regulons | Multi-omics integration | Requires paired RNA-seq; chromVAR alone is insufficient; decoupleR on chromVAR z-scores is not TF activity (motif z-scores are not gene targets) |

## chromVAR versus footprinting

| Question | Tool |
|----------|------|
| Does the bulk pattern of motif-containing peaks vary with condition? | chromVAR |
| Is THIS specific motif site bound by a TF? | TOBIAS / HINT-ATAC |
| Per-cell TF activity in scATAC | chromVAR (via Signac/ArchR) |
| Per-cell TF binding at specific sites | scprinter |
| TF activity correlated with expression | chromVAR + co-expression, or SCENIC+ |
| Which TF families distinguish cell clusters? | chromVAR per-cluster z-scores |
| Differential bound vs unbound between conditions | TOBIAS BINDetect |

## Per-tool failure modes

**chromVAR, too few peaks or reads.** Trigger: peakset < 5000 peaks, or per-sample total depth < 1500 reads. Sparse background sampling creates correlated nulls that inflate positive and negative z-scores. Symptom: variability all > 5; top motifs dominated by AT- or GC-rich sequences regardless of biology. Fix: use a full-scale peakset (typically 50k-200k); for scATAC aggregate cells to clusters of at least 500 cells.

**Custom background peaks.** Trigger: non-default `niterations` or `bias` in `getBackgroundPeaks()`. Default `niterations=50` yields 50 matched background peaks per foreground peak; reducing it adds noise (variability inflated below 30), increasing it slows linearly with little accuracy gain. Fix: keep defaults unless benchmarking; test on a subsample first for huge cell counts.

**Broadly accessible or heterogeneous cell types.** Trigger: cell types with very different overall accessibility. Correction normalizes for total accessibility, so high-background cell types get compressed z-scores. Symptom: PCA on z-scores separates cell types less cleanly than raw counts. Fix: run chromVAR per cell-type cluster when global accessibility differs by more than 5x, or use ArchR's per-cluster background.

**Bulk samples without enough variation.** Trigger: technical replicates or near-identical samples. Z-scores normalize across the population, so with no variability they collapse to zero and variability rankings are unstable. Fix: chromVAR needs about 6+ samples with biological variation; otherwise use footprinting or differential accessibility.

**Signac chromVAR, motif mismatch.** Trigger: motif assay added before the peakset was finalized. Motifs are matched to peaks at call time, so later peak changes leave stale annotations (NA motif annotations, missing deviation entries). Fix: run `AddMotifs()` then the chromVAR step after finalizing peaks; rerun if peaks change.

**ArchR, TileMatrix versus PeakMatrix.** Peaks are biologically meaningful; tiles add intergenic background noise. Use `matrixName='PeakMatrix'` after `addReproduciblePeakSet`; tile deviations are mainly for embedding.

## Reconciling chromVAR, ArchR, and Signac

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Top variable motifs disagree | Different motif databases (JASPAR vs CIS-BP) | Re-run with a matched motif set |
| Z-scores correlate but magnitudes differ | Different background sampling | Inspect per-tool background; defaults are similar but not identical |
| Signac chromvar assay has NA values | Motifs added after peakset finalized | Re-run AddMotifs + chromVAR after peaks are stable |
| ArchR per-cluster signature differs from Signac | Different clustering or cell membership | Standardize clustering before comparison |

Operational rule: chromVAR z-scores are tool-specific. For cross-study comparison, recompute on the same peakset with the same motif database; do not use stored z-scores from heterogeneous sources directly.

## Variability interpretation

| Variability | Typical z-score range | Interpretation |
|-------------|----------------------|----------------|
| < 1 | -1 to +1 | Roughly constant; not biologically variable |
| 1-2 | -2 to +2 | Modest variation; condition-driven possible |
| 2-5 | -3 to +5 | Strong cross-sample / cross-cluster variability; biologically interesting |
| > 5 | -5 to +10 | Major driver of cell-state differences; flagship hits |

Variability is the across-sample SD of z-scores (checked against the bulk run); the ranges above are inherited source guidance, not measured, and depend on the contrast; it ranks motifs without condition labels, so it is the primary metric for unsupervised TF discovery (for example trajectory analysis).

## Background peak matching

chromVAR matches each foreground peak to background peaks by GC content and total accessibility using bin size `bs` (default 50), sampling `niterations` (default 50) replacements from the matching bins; the variance across these samples is the null reference.

- `bs=50` (default) suits typical peaksets; for very small peaksets (< 2000 peaks) lower `bs` to avoid empty bins.
- `niterations=50` (default): below 30 inflates noise; above 100 gives diminishing returns.
- For non-canonical genomes (for example mm10 with different GC distribution), consider rebuilding bins manually with `quantile()` to ensure equal-sized bins.

## chromVAR versus sequence-model alternatives (single cell)

| Tool | Approach | Best for | Limitation |
|------|---------|----------|------------|
| chromVAR | Matched-background z-score per motif | Standard sc workflow; integrated in Signac/ArchR | Linear; no sequence context beyond motif PWM |
| scBasset (Yuan & Kelley 2022) | Sequence CNN with per-cell projection | Higher cluster-discrimination accuracy than chromVAR | Newer; smaller ecosystem; needs >= 100 cells per cluster |
| Enformer-derived TF activity | Long-context Transformer | Cross-cell-type TF activity; distal regulation | Pre-trained models are cell-type-specific |

For high-stakes per-cell TF activity, run chromVAR plus scBasset and report the intersection. See atac-seq/deep-learning-atac for scBasset details.

## Common errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `Error in addGCBias`: missing `seqlengths` | GRanges lacks chrom sizes | `seqlengths(peaks) <- seqlengths(genome)` first |
| All z-scores near zero | Too few samples or too little variation | chromVAR requires biological variation; use footprinting or differential instead |
| `getBackgroundPeaks` slow | Default niterations and large peakset | Default is fine; do not reduce iterations below 30 |
| Differential motifs all significant | FDR not applied, or identical samples compared | Apply BH correction; verify groups |
| Signac chromvar assay all zero | chromVAR run before peakset was final | Re-run after AddMotifs and peakset stability |
| FindAllMarkers reports `avg_log2FC` for chromvar | Default fc method inappropriate for z-scores | Use `mean.fxn=rowMeans` and `fc.name='avg_diff'` |
| z-score interpretation flipped | Contrast sign reversed | Verify factor level order; first level is the reference |
| ArchR `cisbp` vs `JASPAR2020` results differ | Different motif databases | Choose one and report it |

## References

- Schep AN et al 2017 Nat Methods 14:975 (chromVAR)
- Granja JM et al 2021 Nat Genet 53:403 (ArchR)
- Stuart T et al 2021 Nat Methods 18:1333 (Signac)
- Aibar S et al 2017 Nat Methods 14:1083 (SCENIC)
- Bravo Gonzalez-Blas C et al 2023 Nat Methods 20:1355 (SCENIC+)
- Castro-Mondragon JA et al 2022 NAR 50:D165 (JASPAR 2022)
- Rauluseviciute I et al 2024 NAR 52:D174 (JASPAR 2024)
- Weirauch MT et al 2014 Cell 158:1431 (CIS-BP)
- Vorontsov IE et al 2024 NAR 52:D154 (HOCOMOCO v12)
- Yuan H & Kelley DR 2022 Nat Methods (scBasset)
