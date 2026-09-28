---
name: bio-data-visualization-dimensionality-reduction-plots
description: Produce and interpret PCA, t-SNE, UMAP, and PHATE plots for high-dimensional omics data with rigor about which method preserves what (variance, local structure, manifold, transitions), hyperparameter sensitivity, and the well-documented limits of 2D embeddings. Covers PCA biplot/scree/loadings, t-SNE PCA initialization (Kobak-Berens 2019), UMAP n_neighbors/min_dist trade-offs, and the Chari-Pachter 2023 critique. Use when visualizing high-dimensional data — bulk PCA, single-cell embeddings, multi-omics integration projections.
tool_type: mixed
primary_tool: scanpy
---

## Version Compatibility

Reference examples tested with: scanpy 1.10+, anndata 0.10+, scikit-learn 1.4+, umap-learn 0.5+, openTSNE 1.0+, phate 1.0+, ggplot2 3.5+, PCAtools 2.16+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Dimensionality-Reduction Plots

**"Make a PCA / UMAP / t-SNE plot"** -> Choose a projection method aligned with what the plot must reveal — variance explained (PCA), local neighborhood structure (t-SNE), manifold approximation with some global structure (UMAP), or continuous transitions (PHATE). Set hyperparameters deliberately. Communicate the projection's limits and refuse to over-interpret 2D distances.

- Python: `sklearn.decomposition.PCA`, `openTSNE`, `umap-learn`, `phate`, `scanpy.tl.umap` / `scanpy.tl.tsne` / `scanpy.tl.pca`
- R: `prcomp`, `PCAtools::pca`, `Seurat::RunPCA` / `RunUMAP` / `RunTSNE`, `phateR`

## The Single Most Important Modern Insight -- 2D Embeddings Distort

Chari & Pachter 2023 *PLOS Comp Biol* 19:e1011288 demonstrated that 2D embeddings of single-cell data lose >95% of the high-dimensional geometry. Local neighborhoods are only partly preserved (often well under half of the 15 nearest neighbors), while distances between distant cells, density estimates, and global topology are not preserved. The "specious art" of single-cell genomics is the practice of reading 2D layout as biology.

Practical consequence: a UMAP plot communicates "these cells are similar locally" and nothing more. Distance between clusters is meaningless. Density of points within a cluster is dominated by the embedding's repulsion parameter, not the underlying biology. A trajectory inferred from "the gap" between two clusters in UMAP space is an artifact unless validated against the high-dimensional data (RNA velocity, diffusion pseudotime, PHATE).

A second foundational paper is Kobak & Berens 2019 *Nat Commun* 10:5416 on t-SNE for single-cell: PCA initialization + early-exaggeration + multi-scale similarity kernels recover more global structure than default t-SNE settings. The same logic applies to UMAP via `init='spectral'` (default) and `min_dist`.

## Algorithmic Taxonomy

| Method | Preserves | Hyperparameters | Strength | Fails when |
|--------|-----------|-----------------|----------|------------|
| PCA | Linear variance (orthogonal, ordered) | n_components, scaling, solver/seed | Interpretable via loadings; reproducible with a fixed solver or seed; variance % per axis | Non-linear manifolds; high-dim data with few effective dims |
| t-SNE (van der Maaten 2008) | Local neighborhoods (Student-t similarity) | perplexity (typ. 30-50), learning_rate, n_iter, init | Crisp cluster separation | Global distances meaningless; cluster sizes deceptive; non-deterministic |
| UMAP (McInnes 2018, Becht 2018) | Manifold local + partial global | n_neighbors (typ. 15-50), min_dist (typ. 0.1-0.5), spread | Faster than t-SNE; better global preservation than default t-SNE; deterministic given seed | Still distorts; n_neighbors small -> shattered; large -> homogenized |
| PHATE (Moon 2019) | Continuous transitions, branching trajectories | k (knn), t (diffusion power) | Best for developmental trajectories; preserves transition geometry | Slower; less canonical for clustering display |
| Diffusion map | Diffusion distance | epsilon, n_components | Theoretically motivated; supports pseudotime | Less visually striking; less commonly used |
| MDS / classical MDS | Global Euclidean distances | n_components, dissimilarity matrix | Honest about distance preservation | Computationally expensive >5000 points |
| Isomap | Geodesic distance on knn graph | n_neighbors, n_components | Captures non-linear manifold | Sensitive to k; less popular than UMAP |
| Force-directed (PAGA, ForceAtlas2) | Graph topology | Layout-specific | Best for connectivity (PAGA cluster graph) | Not for dense cells; aesthetic |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Bulk RNA-seq sample QC | PCA on log-vst counts; show PC1 vs PC2 with metadata color | Variance explained is meaningful for batch detection |
| Single-cell broad cluster overview | UMAP `n_neighbors=30, min_dist=0.3` after PCA(50) | Standard; preserves clusters; faster than t-SNE |
| Single-cell with delicate trajectories | PHATE OR diffusion map | Preserves continuous transitions |
| Cluster cardinality / boundary visualization | t-SNE with PCA init, perplexity=50 (Kobak-Berens) | Crisper cluster separation than UMAP |
| Multi-omics integration projection | MOFA factors + PCA, or UMAP of joint embedding | Per-omics projection often misleading |
| Spatial transcriptomics with histology | UMAP for transcriptional axis; SEPARATE spatial scatter | UMAP collapses physical space |
| Identify which genes drive variation | PCA biplot with loadings as arrows | Loadings are interpretable; UMAP/t-SNE has no loadings |
| Demonstrating batch confound | PCA colored by batch; quantify batch association across PC1-PC5 | UMAP is a display, not a batch-effect statistic; use PC-wise R2 or kNN batch mixing |
| Visualizing 50 conditions | UMAP/t-SNE for nuance; faceted PCA for interpretability | Method choice depends on question |

## Reference Files

- [Method recipes](references/method-recipes.md): read when implementing PCA/loadings, t-SNE, UMAP, Scanpy saving, or PHATE. The decision-tree rows above map directly to those recipes.
- [Neighborhood validation](references/neighborhood-validation.md): read before claiming that an embedding preserved neighborhoods, batch mixing, or global structure.
- [`examples/embedding_phd.py`](examples/embedding_phd.py): runnable Scanpy comparison with raw-count HVG selection, categorical legends, PCA loading arrows, five saved figures, fixed seeds, and a neighborhood-retention measurement.

## Method Rules

- **PCA:** normalize library size before log/scale (or use `vst`/`rlog`), use a fixed full solver or seed, label axes with variance explained, and inspect loadings. A low PC1/PC2 percentage is a warning for bulk sample QC; it is normal in sparse single-cell data and is not evidence of noise by itself.
- **t-SNE:** use PCA initialization and a learning rate scaled to sample count. In openTSNE 1.0 these are already the `initialization='pca'` and `learning_rate='auto'` defaults; spell them out for provenance, not because openTSNE defaults to random/200. Rtsne has no `seed` argument: call `set.seed(42)` before `Rtsne()`.
- **UMAP:** `min_dist` controls rendered tightness and `n_neighbors` controls the local/global trade-off. Set `random_state` for umap-learn or `seed` for uwot; Scanpy already defaults `random_state=0`, but pass it explicitly for recorded provenance.
- **Scanpy output:** call `sc.pl.umap(..., show=False)` and save the returned axis's figure with `fig.savefig(path, dpi=300)` for exact path control. The older `save=` convention prefixes filenames under `sc.settings.figdir` and is deprecated in Scanpy 1.12.
- **PHATE:** use for continuous or branching trajectories, then validate the inferred ordering against pseudotime, velocity, or known markers.

## Per-Method Failure Modes

### Over-interpreting UMAP distances

**Trigger:** Reading "cluster A is closer to cluster B than to C" as biological similarity.

**Mechanism:** UMAP optimizes a local-neighborhood objective but preserves only part of the original neighbor set; global distances are not preserved (Chari-Pachter 2023).

**Symptom:** Conclusion contradicts hierarchical clustering / RNA velocity / known biology.

**Fix:** Validate inter-cluster relationships against high-dimensional metrics (correlation, distance in PCA space, RNA velocity).

### t-SNE / UMAP without random seed

**Trigger:** Reproducibility request; figure differs between runs.

**Mechanism:** umap-learn, uwot, openTSNE, and Rtsne use stochastic optimization. Scanpy currently supplies a deterministic default (`random_state=0`), but explicit seeds make provenance portable across implementations.

**Symptom:** Re-running the script produces visibly different layouts.

**Fix:** Set `random_state=42` (umap-learn/openTSNE), `seed=42` (uwot), or call `set.seed(42)` immediately before `Rtsne()`; Rtsne silently ignores `seed=` through `...`.

### Perplexity too high for sample size

**Trigger:** Default t-SNE perplexity (30) on small dataset (<500 points).

**Mechanism:** Perplexity >= n/3 is invalid. Rtsne errors; openTSNE clamps it with a warning.

**Symptom:** Plot shows "shattered" small clusters that don't correspond to biology.

**Fix:** For small n: perplexity = max(5, n/30). For very large n: perplexity 50-100.

### PCA without scaling

**Trigger:** `prcomp(X)` or `PCA().fit(X)` without scaling rows/columns first.

**Mechanism:** Genes with high absolute expression dominate variance; PCA captures library-size effect rather than biological variation.

**Symptom:** PC1 perfectly correlates with library size or with total expression.

**Fix:** Normalize for library size before log + scale, or use `vst()` / `rlog()` for bulk counts. Scaling alone does not remove a library-size component. For single-cell, normalize totals, log-transform, then scale.

### UMAP cluster shapes "interpreted" as biological signal

**Trigger:** Reporting that a cluster is "elongated" or "round" as biological observation.

**Mechanism:** UMAP cluster shape is an artifact of `min_dist` and `n_neighbors`, not biology.

**Symptom:** Reviewer asks "why is the immune cluster elongated?"; no answer except the embedding.

**Fix:** Do not interpret cluster shape. Report cluster membership and validate biology via marker genes.

### scanpy.pl.umap save writes to figures/ subdirectory

**Trigger:** `sc.pl.umap(adata, save='myplot.pdf')` with the expectation that myplot.pdf will land in the current directory.

**Mechanism:** `save=` is concatenated with `sc.settings.figdir` (default `./figures/`) and prefixed with `umap`.

**Symptom:** File not at the requested path; actually at `figures/umapmyplot.pdf`.

**Fix:** Set `sc.settings.figdir='/abs/path/'` AND `save='_descriptive.pdf'` so result is `figures/umap_descriptive.pdf`. For full path control use `matplotlib.savefig` after `sc.pl.umap(show=False)`.

### Default DPI 150 below journal requirements

**Trigger:** scanpy `sc.settings.set_figure_params()` default `dpi_save=150`.

**Mechanism:** Nature/Cell require 300+ DPI for raster.

**Symptom:** Figure looks fine on screen, rejected at submission.

**Fix:** `sc.set_figure_params(dpi_save=300, figsize=(4, 4))`.

### PCA loadings interpreted on UMAP/t-SNE coordinates

**Trigger:** "PC1 axis on UMAP" — projecting loadings onto UMAP.

**Mechanism:** UMAP/t-SNE coordinates have no linear interpretation; loadings are PCA-specific.

**Symptom:** Conclusion about "what UMAP-x means" that has no foundation.

**Fix:** Use PCA when loadings are needed. Show UMAP for visualization and PCA for axis-driving gene identification, separately.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| t-SNE shows distinct clusters; UMAP merges them | t-SNE over-emphasizes local structure; UMAP n_neighbors too large | Both views valid; check Leiden cluster assignments rather than embedding |
| A top PC is associated with batch; UMAP looks mixed or separated | Embedding appearance does not quantify batch | Screen PC1-PC5 (or more) by batch R2 and measure kNN batch mixing; correct upstream before biological interpretation |
| PHATE shows continuous trajectory; UMAP shows discrete clusters | PHATE preserves transitions; UMAP "blobifies" continuous data | Use PHATE for trajectory display; UMAP for discrete cell-type display |
| Reproducibility breaks across re-runs | Random seed not set | Set seed; document version of umap-learn/openTSNE |
| Cluster boundaries differ between Seurat/scanpy UMAP | Different defaults for n_neighbors, min_dist, init | Standardize hyperparameters; report explicitly |

**Operational rule:** report ALL hyperparameters used (perplexity, n_neighbors, min_dist, random_state). State the embedding's interpretation limit ("local neighborhood; distances between clusters not meaningful"). For trajectory claims, validate with RNA velocity, pseudotime, or PHATE.

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| t-SNE perplexity | 30-50 for n>1000; max(5, n/30) for n<500 | Maaten 2008; Kobak-Berens 2019 |
| UMAP n_neighbors | 15-50 default range | umap-learn docs; Becht 2018 |
| UMAP min_dist | 0.1-0.5 | Tighter for crisp clusters, looser for continua |
| t-SNE learning rate | n / 12 (Kobak-Berens; openTSNE `auto`) | Set explicitly outside openTSNE when defaults differ |
| PCA n_components for downstream UMAP | 30-50 | Standard scanpy workflow |
| Single-cell n_neighbors for sc.pp.neighbors | 15-30 | Wolf 2018 Scanpy paper |
| Save DPI for publication | 300+ | Nature/Cell figure guidelines |
| Random seed | 42 (or any fixed integer) | Reproducibility |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Plot differs between runs | Random seed not fixed in the chosen implementation | Record an explicit supported seed; use `set.seed()` for Rtsne |
| PC1 = library size | No library-size normalization | Normalize totals before log + scale, or use `vst()`/`rlog()` |
| Cluster shapes "interpreted" biologically | UMAP artifact | Do not interpret shape; report membership |
| scanpy save writes to wrong path | figdir + prefix concatenation | Set figdir explicitly OR use matplotlib.savefig |
| t-SNE fragments small dataset | Perplexity too high for n | Use perplexity = max(5, n/30) |
| Inter-cluster "distance" used in trajectory claim | UMAP distances not meaningful | Validate with RNA velocity / PHATE |
| Loadings interpreted on UMAP axes | UMAP has no loadings | Use PCA for loading-driven interpretation |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why UMAP and not t-SNE?" | UMAP for cluster overview (faster, better global preservation given Becht 2018); t-SNE in supplementary if cluster boundaries are the focus |
| "What hyperparameters?" | Explicit n_neighbors, min_dist, n_pcs, random_state in caption AND methods |
| "Why is cluster X shaped this way?" | UMAP/t-SNE cluster shape is an embedding artifact; cluster membership is the biological observation |
| "Are these trajectories real?" | Validated via RNA velocity / PHATE / diffusion pseudotime (NOT inferred from UMAP layout alone) |
| "Why PCA?" | Variance explained per axis is interpretable for sample QC; loadings identify driving genes (uniquely PCA, NOT UMAP) |

## References

- Becht E, McInnes L, Healy J, et al. 2019. Dimensionality reduction for visualizing single-cell data using UMAP. *Nat Biotechnol* 37(1):38-44.
- Chari T, Pachter L. 2023. The specious art of single-cell genomics. *PLOS Comp Biol* 19(8):e1011288.
- Kobak D, Berens P. 2019. The art of using t-SNE for single-cell transcriptomics. *Nat Commun* 10:5416.
- McInnes L, Healy J, Melville J. 2018. UMAP: Uniform Manifold Approximation and Projection for dimension reduction. *arXiv:1802.03426*.
- Moon KR, van Dijk D, Wang Z, et al. 2019. Visualizing structure and transitions in high-dimensional biological data. *Nat Biotechnol* 37(12):1482-1492.
- van der Maaten L, Hinton G. 2008. Visualizing data using t-SNE. *J Mach Learn Res* 9:2579-2605.
- Wolf FA, Angerer P, Theis FJ. 2018. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* 19:15.

## Related Skills

- single-cell/preprocessing - PCA / neighbor-graph computation before embedding
- single-cell/clustering - Leiden / Louvain cluster assignments visualized in UMAP
- single-cell/trajectory-inference - PHATE, diffusion pseudotime, RNA velocity for trajectory claims
- data-visualization/color-palettes - Categorical palette for cluster labels
- data-visualization/distribution-plots - Per-cluster gene-expression follow-up
- data-visualization/heatmaps-clustering - Alternative view of the same high-dim matrix
