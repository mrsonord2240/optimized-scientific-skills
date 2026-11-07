# Normalization, HVGs, scaling

Defaults in the scripts: shifted-log, then `seurat_v3` HVGs on raw counts, then PCA on log-normalized HVGs with no scaling. Change only if the request says to.

```python
adata.layers['counts'] = adata.X.copy()   # stash raw before normalizing
sc.pp.normalize_total(adata)              # target_sum=None scales to the dataset median; 1e4 (CP10k) is the arbitrary legacy choice
sc.pp.log1p(adata)
```

| Normalization | Use when | Fails when |
|---------------|----------|-----------|
| Shifted-log (median or CP10k) | General default | Composition-divergent types (plasma cells, cycling cells with 2-10x the mRNA) distort fold changes |
| scran deconvolution (R) | Low depth, high dropout, plate-based | Needs `quickCluster`; size factors can go negative |
| sctransform v2 (Seurat) | Depth removal for HVGs and visualization | Slow; off the count scale |
| Analytic Pearson residuals | UMI data, HVGs plus PCA without ad hoc steps | Experimental; depends on theta and depth; clip to +/- sqrt(n) |

Plain shifted-log plus PCA performed as well as or better than sctransform, Pearson residuals and GLM-PCA on kNN-overlap recovery, so the simple method is the defensible default. Single-cell measures proportions: report relative, not absolute, expression.

| HVG flavor | Function | Input | `n_top_genes` |
|------------|----------|-------|---------------|
| `seurat` (default) | `sc.pp.highly_variable_genes` | log-normalized | optional |
| `cell_ranger` | same | log-normalized | optional |
| `seurat_v3` | same | RAW counts (`layer='counts'`) | required; needs scikit-misc |
| `pearson_residuals` | `sc.experimental.pp.highly_variable_genes` | RAW counts | recommended |

- Set `batch_key` to compute HVGs per batch. With `seurat_v3` it also needs genes detected in at least `min_cells` cells in every batch, or a shallow batch fails with `ValueError: reciprocal condition number`. `scripts/preprocess_scanpy.py batch_key=<col>` applies that gene filter.
- Scaling and regressing out: `sc.pp.scale(adata, max_value=10)` is optional and increasingly skipped (it inflates noisy low-expression genes). `regress_out` of `total_counts` or `pct_counts_mt` erases biology and Scanpy warns it overcorrects. Address unwanted variation with integration (Harmony, scVI).
- Exposing a dominant gene spike: `exclude_highly_expressed=True` in `normalize_total`, or scran.

Done when the normalization and HVG flavor used are stated with their inputs.
