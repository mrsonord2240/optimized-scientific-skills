# Methods and interpretation

Source performance and API claims were not independently verified. A physical multiplexed capture can contain cells from different samples; use the capture as the grouping unit and total loading for rate assumptions.

## Expected Doublet Rate

10X Chromium loading is near-Poisson, so the multiplet rate scales roughly linearly with recovered cells: ~0.8% per 1,000 cells recovered (`dbr.per1k = 0.008`).

| Cells recovered (~) | Expected rate (~) |
|---------------------|-------------------|
| 1,000 | 0.8% |
| 2,000 | 1.6% |
| 5,000 | 3.9% |
| 10,000 | 7.6-8% |

Rule: rate ~= 0.008 x recovered/1000. Always set the expected rate from the actual recovered-cell count of that lane; Scrublet's flat `expected_doublet_rate=0.05` is a placeholder, not a recommendation. High-throughput chips have lower per-cell rates. For multiplexed pools (genotype/HTO-demultiplexed), the physical doublet rate is set by TOTAL lane loading, not the demultiplexed subset: deriving the rate from one sample's cells underestimates it (four 5k samples in one 20k lane is ~15% real, not the ~3.9% implied by 5k), so set the rate from the total lane cell count.

## Heterotypic, Homotypic, Neotypic

| Type | Composition | Detectability |
|------|-------------|---------------|
| Heterotypic | Two transcriptionally distinct types | Detectable; lands between clusters; the dangerous "fake transitional" ones |
| Homotypic | Two cells of the same type | Nearly undetectable by expression; persists after removal |
| Neotypic | Heterotypic blend occupying a region no singlet occupies | Most detectable; most misleading if missed (looks like a rare new type) |

`modelHomotypic`-style adjustments only change the number expected to be detectable; they cannot recover undetectable homotypic doublets.

## Choosing a Method

| Method | Model | Use when | Fails / weak when |
|--------|-------|----------|-------------------|
| scDblFinder (R) | xgboost on kNN features vs simulated doublets | Default; best accuracy-speed balance; built-in per-sample via `samples=` | R/Bioconductor only |
| Scrublet (Python) | kNN density of simulated doublets | scanpy-native pipelines | Auto-threshold fails on unimodal histograms; loop per sample manually |
| DoubletFinder (R) | pANN from PC neighborhood | Legacy Seurat workflows | Brittle; `pK` needs per-dataset sweep; `*_v3` names removed; Seurat-version-coupled |
| solo (Python) | scVI VAE + classifier | Have a trained scVI model; GPU available | Heavier setup |
| scds cxds/bcds/hybrid (R) | Co-expression / boosted tree | Fast first pass on very large data | Lower accuracy than scDblFinder/DoubletFinder |

scDblFinder is the 2024-2026 best-balance default and is recommended by sc-best-practices. The older Xi and Li 2021 ranking ("DoubletFinder is most accurate") predates major scDblFinder improvements and is superseded - do not cite it against current scDblFinder. Methods compete and drift; verify current standing against the installed tool's docs before committing.

## DoubletFinder (R, legacy Seurat)

**Goal:** Run DoubletFinder on a fully preprocessed Seurat object, tuning `pK` and homotypic-adjusting the expected count.

**Approach:** Sweep `pK`, pick the BCmvn maximum, then set `nExp` from the rate adjusted for the homotypic fraction.

```r
library(DoubletFinder)                                     # *_v3 function names were removed in Nov 2023; verify installed API

sweep.res <- paramSweep(seurat_obj, PCs = 1:20, sct = FALSE)
bcmvn <- find.pK(summarizeSweep(sweep.res, GT = FALSE))
pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))   # no default pK; tune per dataset

rate <- 0.008 * ncol(seurat_obj) / 1000
nExp <- round(rate * ncol(seurat_obj))
nExp <- round(nExp * (1 - modelHomotypic(seurat_obj$seurat_clusters)))  # discount undetectable homotypic doublets
seurat_obj <- doubletFinder(seurat_obj, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp, sct = FALSE) # omit reuse.pANN unless supplying an existing pANN column name
```

`pN` (artificial-doublet proportion) defaults to 0.25 and performance is largely insensitive to it. DoubletFinder requires a normalized, PCA'd, clustered object and is the most version-sensitive of the three.

## Deeper Cautions

Simulated doublets are a model, not the real thing: real doublets share one RT/PCR reaction (capture competition, barcode effects), so simulated-doublet density only approximates where real doublets sit, and even the best method has a low ceiling (max mean AUPRC ~0.537 in Xi and Li 2021 - every method misses a lot). Over-removal culls proliferating (S/G2M) and genuine transitional cells that legitimately score high, so cross-check removed cells against cell-cycle and activation signatures. When available, experimental ground truth beats inference: cell hashing (CITE-seq HTOs) and MULTI-seq call inter-sample doublets directly regardless of expression similarity (catching even cross-sample homotypic doublets), and serve as a complementary filter. Heavy ambient RNA can mimic co-expression and nudge scores, so handle empty droplets and ambient RNA first (see single-cell/preprocessing).

## What to report

- Per sample/capture: cells scored, expected doublet rate, and whether that rate came from its own recovered count or the total multiplexed lane.
- Method, tool version, seed/BPPARAM (when applicable), and whether the threshold was automatic or manually chosen after histogram review.
- Doublet calls as both a count and percentage, retained as flags/scores unless removal is explicitly requested.
- The homotypic limitation: reported calls estimate the expression-detectable fraction, not a doublet-free population.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Doublet calls look random / too many | Run on merged multi-sample data | Run per sample before integration (`samples=` or loop) |
| `doublet_score` / `predicted_doublet` are missing after a per-sample loop | Scrublet was run on an AnnData view, so results were not written to the parent | Run each subset with `.copy()`, assign both columns back with `.obs.loc[...]`, and stop if any scores remain missing |
| Auto-threshold splits the histogram badly | Scrublet histogram is unimodal | Inspect the histogram and set `threshold` manually |
| A high-RNA cell type was wiped out | Count-based QC and doublet removal double-penalized the same axis | Coordinate the filters; do not stack aggressive cutoffs |
| "Novel transitional state" co-expresses two lineages | Heterotypic doublets masquerading as a cluster | Confirm per-sample detection; check marker co-expression / hashing before claiming a new type |
| Trajectory has an implausible bridge between lineages | Doublets forming a false intermediate | Remove/flag doublets before trajectory inference |
| Reported "0% doublets" | Homotypic doublets are invisible | Do not claim doublet-free; report only the detectable fraction |
| DoubletFinder errors with `cannot xtfrm data frames` | `reuse.pANN = FALSE` is invalid: the argument is a prior pANN column name or omitted | Omit `reuse.pANN` unless reusing a named pANN column; then re-tune `pK` |
| Expected rate clearly wrong | Used a package default | Set rate from recovered cells (~0.008 x cells/1000) |
| Multiplexed pool underestimates doublets | Rate derived from one demultiplexed sample, not total lane | Set the expected rate from total capture-lane cells |

## References

- Wolock SL, Lopez R, Klein AM (2019) Scrublet: computational identification of cell doublets in single-cell transcriptomic data. Cell Systems 8(4):281-291.e9. DOI 10.1016/j.cels.2018.11.005
- McGinnis CS, Murrow LM, Gartner ZJ (2019) DoubletFinder: doublet detection in single-cell RNA sequencing data using artificial nearest neighbors. Cell Systems 8(4):329-337.e4. DOI 10.1016/j.cels.2019.03.003
- Germain P-L, Lun A, Macnair W, Robinson MD (2021) Doublet identification in single-cell sequencing data using scDblFinder. F1000Research 10:979. DOI 10.12688/f1000research.73600
- Xi NM, Li JJ (2021) Benchmarking computational doublet-detection methods for single-cell RNA sequencing data. Cell Systems 12(2):176-194.e6. DOI 10.1016/j.cels.2020.11.008
- Bernstein NJ, Fong NL, Lam I, et al. (2020) Solo: doublet identification in single-cell RNA-seq via semi-supervised deep learning. Cell Systems 11(1):95-101.e5. DOI 10.1016/j.cels.2020.05.010
- Bais AS, Kostka D (2020) scds: computational annotation of doublets in single-cell RNA sequencing data. Bioinformatics 36(4):1150-1158. DOI 10.1093/bioinformatics/btz698
- McGinnis CS, Patterson DM, Winkler J, et al. (2019) MULTI-seq: sample multiplexing for single-cell RNA sequencing using lipid-tagged indices. Nature Methods 16(7):619-626. DOI 10.1038/s41592-019-0433-8
- Heumos L, Schaar AC, Lance C, et al. (2023) Best practices for single-cell analysis across modalities. Nature Reviews Genetics 24:550-572. DOI 10.1038/s41576-023-00586-w
