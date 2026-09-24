---
name: bio-single-cell-batch-integration
category: Data Analysis
description: Integrate multiple scRNA-seq samples or batches with Harmony, scVI/scANVI, Seurat (CCA/RPCA), fastMNN, Scanorama, or BBKNN. Resolves which method to use for the dataset size and design, how strongly to correct, when integration is the wrong move (confounded batch/biology), how to score integration with scIB metrics without gaming them, and why corrected expression must not be used for differential expression. Use when integrating batches or datasets, choosing an integration method, diagnosing over-correction, or judging integration quality.
tool_type: mixed
primary_tool: Harmony
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: scanpy 1.10+, Seurat 5.0+, scvi-tools 1.1+, harmonypy 0.0.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Batch Integration

**"Integrate my batches"** -> Learn a shared low-dimensional representation that mixes technical batches while preserving biological cell states, then cluster and visualize on it.
- Python: `sce.pp.harmony_integrate`, `scvi.model.SCVI`, `sce.pp.bbknn`, `scanorama`
- R: `RunHarmony`, `IntegrateLayers` (Seurat v5), `fastMNN` (batchelor)

## Governing Principle

Integration trades batch-mixing against biological-signal preservation, and the two cannot be jointly maximized. The algorithm removes variance along directions where batches differ, on the assumption that cell-type composition is shared across batches; it does not "know" which variance is technical. When batch correlates with a real biological axis, the method cannot distinguish them and, by construction, erases biology - this is information-theoretic, not a tuning problem. Over-correction is the silent failure: it absorbs rare cell types into common neighbors, snaps continuous gradients toward shared anchors, and deletes condition-specific populations, all while batch-mixing metrics improve. The cells most worth finding (rare, transitional, novel) are exactly the ones integration most endangers. Batch and biology are unidentifiable under confounding - the only tie-breaker is external information (shared controls, multiplexed designs, known shared types), and the fix for a confounded design is experimental, not computational. Always keep the uncorrected embedding for before/after comparison, and never run differential expression on batch-corrected expression.

## When NOT to Integrate

Visualize the uncorrected data first; integration is a bias-variance trade and removing batch variance risks removing biology correlated with batch.

- Confounded design (each condition is its own batch, e.g. all controls day 1, all treated day 2): no algorithm can separate batch from biology. **Before looking at an embedding, cross-tabulate condition x batch:** stop if any condition occurs in only one batch. Then cluster the uncorrected data and cross-tabulate clusters x batch x condition to size a visible technical effect; this follow-up does not make a confounded design safe when clusters mix. The fix is experimental - multiplex conditions across batches (cell hashing; genetic demux via souporcell/vireo; split each condition across capture days).
- Technical replicates of the same tissue that already mix well: over-correction risk outweighs benefit.
- Per-sample analyses (CNV/tumor-clone inference): integration would erase the signal of interest.

Over-correction signatures: rare types collapsing into neighbors, lost known gradients, disappearing condition-specific populations, markers no longer separating known cell types.

## Method Selection

No method wins universally - the scIB benchmark (Luecken 2022, 68 method/preprocessing combos) scores integration as overall = 0.6 x bio-conservation + 0.4 x batch-removal, deliberately weighting biology higher because erasing it is worse than imperfect mixing. Methodology evolves; verify current best practice and APIs against the installed package docs, and in practice run 2-3 candidates and score them (see Evaluating Integration).

| Method | Model / assumption | Use when | Fails when |
|--------|-------------------|----------|------------|
| [Harmony](references/harmony.md) | Iterative soft k-means linear correction in PCA space; outputs an embedding, not counts | Few/simple batches, fast, low memory; strong default; best usability | Strong nonlinear batch effects; high `theta` over-mixes and collapses distinct types |
| [scVI](references/scvi-scanvi.md) | Conditional VAE on raw counts (ZINB), batch as covariate -> batch-invariant latent | Large atlases, many nested batches, strong effects; memory-efficient at scale | Small data (under-trained); latent dims over-interpreted as "biology minus batch" |
| [scANVI](references/scvi-scanvi.md) | Semi-supervised scVI using partial labels to protect biology | Some cell labels exist and bio fidelity is paramount (tops bio-conservation) | Labels noisy/wrong; training cost; closed-world for the labeled states |
| [Seurat CCA](references/seurat-v5.md) | Anchor-based, canonical correlation across datasets | Strong shared structure under large shifts; smaller data | Substantial non-overlap or many samples -> over-correction (CCA aligns distinct states) |
| [Seurat RPCA](references/seurat-v5.md) | Reciprocal-PCA anchors; faster, more conservative | Large/many-sample data, substantial non-overlap | Under-correction when truly shared structure is subtle (raise `k.anchor`) |
| fastMNN | Mutual nearest neighbors in PCA space | Rare-population preservation; moderate data | Order-sensitive (set `merge.order`, most-heterogeneous first); legacy mnnCorrect is slow |
| Scanorama | Mutual NN across all dataset pairs | Partial cell-type overlap across datasets; balanced bio/batch; sort cells by batch key before calling `scanorama_integrate` | Very large data (slower than Harmony/BBKNN) |
| BBKNN | Modifies only the neighbor graph (batch-balanced kNN) | Speed; only clustering/UMAP needed downstream | Leans toward batch removal; no embedding or corrected counts for other uses |

scIB headline: top combined performers were scANVI, scVI, Scanorama, scGen; Harmony and Seurat were strong on simpler tasks with the best usability; BBKNN sits at the batch-removal end. "Deep methods are always best" is not supported - Harmony/Seurat win simple/small tasks; deep methods win complex/large/label-rich tasks.

## Strength Parameters

Aggressive settings increase mixing and over-correction risk in lockstep - raise correction strength only after confirming under-correction, and re-check rare populations after each change.

| Parameter | Tool | Effect | Rationale |
|-----------|------|--------|-----------|
| theta | Harmony | Higher -> more aggressive batch mixing | Start near 2 (roughly 0.5-5 is usually the useful tuning range); below that, batch-effect size can dominate the setting. `theta = 0` does not disable correction, so judge each change with paired batch/bio metrics and rare-population checks. |
| k.anchor | Seurat | Higher -> more anchors, stronger correction | Raise (e.g. 20) only when under-correcting |
| CCA vs RPCA | Seurat | CCA more sensitive but can over-correct; RPCA conservative | Prefer RPCA for large/non-overlapping data |
| n_latent | scVI | Latent dimensionality of the embedding | ~10-30; too high refits noise, too low under-fits |
| merge.order | fastMNN | Order batches are merged | Order-sensitive; merge most-heterogeneous batch first |

## Reference Files

- [`references/harmony.md`](references/harmony.md): read when using Harmony or adjusting `theta`.
- [`references/scvi-scanvi.md`](references/scvi-scanvi.md): read when using scVI/scANVI and preparing raw-count layers or labels.
- [`references/seurat-v5.md`](references/seurat-v5.md): read when integrating a Seurat v5 object, including the RPCA memory setting and optional wrapper methods.

## Integrate with Harmony

**Goal:** Correct batch in PCA space and run downstream steps on the corrected embedding.
**Approach:** Joint preprocessing -> PCA -> Harmony -> neighbors/UMAP/clustering on `X_pca_harmony`.

Run `python scripts/harmony_integration.py merged.h5ad integrated_harmony.h5ad --batch-key batch --seed 0`.
It preserves `X_pca` for the required before/after comparison and writes `X_pca_harmony` for downstream work.

In Seurat: `RunHarmony(obj, group.by.vars = 'orig.ident', reduction.use = 'pca')` writes a `harmony` reduction; `group.by.vars` takes a vector to correct multiple covariates.

## Integrate with scVI / scANVI

**Goal:** Learn a batch-invariant latent space from raw counts, optionally protecting known labels.
**Approach:** Put raw counts in a layer, register batch (and labels for scANVI), train, and use the latent embedding downstream.

Run `python scripts/scvi_scanvi_integration.py merged.h5ad integrated_scvi.h5ad --batch-key batch --labels-key cell_type --max-epochs 100 --seed 0`.
The script keeps raw counts in `layers['counts']`, sets `scvi.settings.seed`, fixes training duration, and writes `X_scVI` (plus `X_scANVI` when `--labels-key` contains `Unknown`).

The scVI latent space is not "biology with batch removed": it is a learned nonlinear embedding optimized to reconstruct counts while being marginally independent of batch. Its dimensions are entangled, individually uninterpretable, and carry no guaranteed correspondence to any biological quantity - treat it as a coordinate system for neighbors/clustering, not a measurement. Note `unlabeled_category` ('Unknown') is the second positional argument to `from_scvi_model`, before `labels_key`.

## Integrate with Seurat v5

**Goal:** Use Seurat v5's modular layer-based integration with a chosen method.
**Approach:** Split layers by batch, run the standard pipeline, call IntegrateLayers, rejoin.

Run `Rscript scripts/seurat_v5_integration.R merged.rds integrated_rpca.rds batch RPCAIntegration` (from this skill directory). The script sets `future.globals.maxSize = 4 * 1024^3` before `IntegrateLayers`; raise it only when the object and available memory justify it.

Seurat supplies `CCAIntegration`, `RPCAIntegration`, and `HarmonyIntegration`. `FastMNNIntegration` and `scVIIntegration` are optional **SeuratWrappers** methods, not bare Seurat symbols; the latter also needs a working reticulate/scvi-tools installation. If SeuratWrappers is unavailable, use `batchelor::fastMNN` or the Python scVI path rather than retrying an unresolved method. For graph-only correction with BBKNN in Python: `sce.pp.bbknn(adata, batch_key='batch')` rewrites the neighbor graph in place (very fast, feeds Leiden/UMAP only).

## Evaluating Integration

**Goal:** Decide whether integration mixed batches without erasing biology.
**Approach:** Score batch-mixing and bio-conservation separately and read them jointly - never optimize a batch metric alone.

```python
import scanpy as sc
from sklearn.metrics import silhouette_score

# batch silhouette: lower = batches mixed; cell-type silhouette: higher = biology kept
batch_sil = silhouette_score(adata.obsm['X_scVI'], adata.obs['batch'])
ct_sil = silhouette_score(adata.obsm['X_scVI'], adata.obs['cell_type'])

# scib-metrics Benchmarker scores many methods on a common axis set
# from scib_metrics.benchmark import Benchmarker
```

Batch-mixing metrics (kBET, graph iLISI) are trivially maximized by over-correction - a method that destroys all structure mixes batches perfectly while annihilating biology. Bio-conservation metrics (ARI, NMI, cell-type ASW, graph cLISI, isolated-label F1) guard against that, which is why the scIB composite down-weights batch-removal to 0.4. Selecting a method on a batch metric alone selects for over-correction; always pair batch metrics with bio metrics and inspect rare populations before/after. Run candidates through `scib-metrics` (`Benchmarker`) and pick the most robust for the specific task.

Differential expression: use integration outputs (Harmony/scVI/RPCA embeddings) for clustering and visualization, but run DE on uncorrected, log-normalized counts - never on batch-corrected expression. Harmony and BBKNN produce no corrected counts; Scanorama and fastMNN do, and those must not feed DE. For cross-condition DE, aggregate to pseudobulk per sample x cell type (see differential-expression/deseq2-basics).

## Reference Mapping vs De-novo Integration

De-novo integration jointly embeds all datasets symmetrically (everything above). Reference mapping projects a query onto a fixed reference embedding without retraining (scArches architectural surgery; Azimuth `FindTransferAnchors` + `MapQuery`) - fast, reproducible, scales to millions, consistent cross-study labels. It is closed-world: a novel state the reference never saw is confidently assigned the nearest reference label, converting a technical or biological surprise into a wrong annotation that looks clean and high-confidence. Use reference mapping when a high-quality annotated atlas exists; use de-novo when no suitable reference exists or the query may hold genuinely novel populations. Always inspect per-cell mapping uncertainty and never trust transferred labels for clusters that map poorly.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| The cell type of interest vanished after integration | Over-correction absorbed a rare/condition-specific population | Reduce strength (lower theta / use RPCA / fastMNN/Scanorama); compare to uncorrected embedding |
| Batches still separate on UMAP | Under-correction | Raise correction strength (k.anchor, switch CCA, more Harmony iterations); confirm batch key is correct |
| "Integration removed my treatment effect" | Confounded batch/condition design | Stop - batch and biology are unidentifiable; redesign with multiplexing; do not integrate away the contrast |
| Great iLISI/kBET but biology looks flattened | Metric gaming by over-correction | Score bio-conservation too (ASW-celltype, cLISI, ARI); use the scIB composite, not a batch metric alone |
| DE between two control samples after integration | DE run on batch-corrected expression | Run DE on uncorrected log-normalized counts; use pseudobulk for cross-condition |
| scVI latent dimension interpreted as a biological axis | Latent space is entangled, not "biology minus batch" | Use the embedding only for neighbors/clustering; do not read individual dims |
| Reference-mapped labels look confident but wrong | Closed-world projection of a novel/shifted state | Inspect mapping uncertainty; treat poorly-mapping clusters as candidate novelty/batch |
| Results differ run to run | Stochastic training / unpinned seeds (scVI, Harmony) | Set seeds; for scVI fix `max_epochs` and report it |
| `future` globals exceed 500 MiB during Seurat RPCA | Anchor integration exports a large object to future workers | Set `options(future.globals.maxSize = 4 * 1024^3)` before `IntegrateLayers`, after checking available memory |
| `Detected non-contiguous batches` from Scanorama | Cells are not contiguous by the batch key | Stable-sort cells by batch before `scanorama_integrate`; retain the original order mapping for plotting/metadata joins |

## Related Skills

- preprocessing - QC and normalization that must precede integration
- clustering - Cluster on the integrated embedding, not on raw PCA
- cell-annotation - Reference mapping and label transfer after integration
- single-cell/multimodal-integration - Joint analysis across modalities (distinct from batch integration)
- single-cell/differential-abundance - Test whether composition shifts across conditions after integration
- single-cell/cnv-inference - Per-patient malignant-cell/CNV inference (do not integrate tumors across patients)
- differential-expression/deseq2-basics - Pseudobulk DE on uncorrected counts per cell type
- data-visualization/dimensionality-reduction-plots - Before/after UMAP comparison figures

## References

- Korsunsky et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods 16(12):1289-1296.
- Lopez et al. (2018). Deep generative modeling for single-cell transcriptomics (scVI). Nat Methods 15(12):1053-1058.
- Xu et al. (2021). Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models (scANVI). Mol Syst Biol 17(1):e9620.
- Luecken et al. (2022). Benchmarking atlas-level data integration in single-cell genomics (scIB). Nat Methods 19:41-50.
- Hie, Bryson & Berger (2019). Efficient integration of heterogeneous single-cell transcriptomes using Scanorama. Nat Biotechnol 37:685-691.
- Polanski et al. (2020). BBKNN: fast batch alignment of single cell transcriptomes. Bioinformatics 36(3):964-965.
- Haghverdi et al. (2018). Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors (MNN). Nat Biotechnol 36:421-427.
- Lotfollahi et al. (2022). Mapping single-cell data to reference atlases by transfer learning (scArches). Nat Biotechnol 40(1):121-130.
