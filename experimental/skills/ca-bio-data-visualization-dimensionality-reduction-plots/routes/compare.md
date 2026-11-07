# Compare embeddings

```bash
python examples/embedding_phd.py raw_counts.h5ad --output-dir figures
```

Input must contain at least 100 cells and 2,000 genes, raw nonnegative integer counts in `X`, and `obs['condition']`; optional `obs['pseudotime']` colors PHATE. The script selects 2,000 highly variable genes before normalizing, then log-transforms and scales; it writes PCA, scree, UMAP, t-SNE and PHATE PDFs and prints a 15-neighbor retention caption. It uses PCA 50, neighbors 30, UMAP min_dist 0.3, t-SNE perplexity 30 and PCA initialization, PHATE k=10, decay=40 and seed 42. It runs all methods, requires their dependencies, and may overwrite named output files; use a new output directory. Category display is limited to 20 and colors are not a substitute for accessibility review. For one method or different input, adapt [method recipes](../references/method-recipes.md); no standalone individual-method CLI is supplied.

Done when plots, parameters, input provenance, category labels, neighbor retention and limitations are reviewed.
