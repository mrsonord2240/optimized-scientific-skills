# Method recipes

Read this file only for the method selected by the decision tree in `SKILL.md`.

## PCA

For bulk count data, normalize library size before PCA; DESeq2 `vst()` or `rlog()` is preferred. PCAtools supplies a biplot, scree plot, and loadings plot:

```r
vsd <- DESeq2::vst(dds, blind = FALSE)
p <- PCAtools::pca(SummarizedExperiment::assay(vsd), metadata = as.data.frame(SummarizedExperiment::colData(dds)))
PCAtools::biplot(p, colby = "condition", shape = "batch", showLoadings = TRUE)
PCAtools::screeplot(p, components = seq_len(min(10, ncol(SummarizedExperiment::assay(vsd)))))
PCAtools::plotloadings(p, components = 1, rangeRetain = 0.05)
```

With scikit-learn, use a fixed solver or seed. Convert string groups to categories and draw one scatter per group so the legend is truthful. The shipped example gives each of up to 20 categories a distinct `tab20` color; above 20, facet the plot or choose and document another encoding instead of recycling colors. `components_` contains feature loadings; the example scales and draws its largest arrows, then deterministically checks rendered label boxes to avoid collisions.

```python
pca = PCA(n_components=10, svd_solver="full")
scores = pca.fit_transform(X_normalized)
loadings = pca.components_.T
```

Label axes from `explained_variance_ratio_`. For batch diagnosis, regress the batch indicator against PC1-PC5 rather than stopping at the first two displayed PCs.

## t-SNE

Kobak-Berens settings improve global organization relative to random initialization in implementations whose defaults differ. openTSNE 1.0 already defaults to PCA initialization and `learning_rate="auto"` (n/12); setting them explicitly records the analysis contract.

```python
embedding = openTSNE.TSNE(perplexity=30, n_iter=750, initialization="pca",
                          learning_rate="auto", random_state=42).fit(X)
```

Rtsne has no `seed` argument. Call `set.seed(42)` before it. For small samples require `perplexity < n/3`; Rtsne errors above that bound while openTSNE warns and clamps.

## UMAP and Scanpy

For umap-learn, record `n_neighbors`, `min_dist`, metric, and `random_state`. For uwot, pass `seed`. Scanpy currently defaults to `random_state=0`, but supply it explicitly so provenance does not depend on an implementation default.

```python
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50)
sc.tl.umap(adata, min_dist=0.3, random_state=42)
ax = sc.pl.umap(adata, color="leiden", show=False)
ax.figure.savefig("umap_clusters.pdf", dpi=300, bbox_inches="tight")
```

## PHATE

Use PHATE for continuous transitions, not as evidence of a trajectory by itself:

```python
emb = phate.PHATE(knn=10, decay=40, t="auto", random_state=42).fit_transform(X)
```

Validate the ordering against pseudotime, velocity, known markers, or experimental time.
