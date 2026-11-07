# Evaluate interpretation

```bash
python examples/embedding_phd.py raw_counts.h5ad --output-dir validation-figures
```

This full comparison prints overlap of the 15 nearest neighbors in PCA and UMAP space; see [comparison](compare.md) for the input contract and overwrite warning. It does not calculate batch association or validate trajectories. Adapt the independently supplied [neighborhood validation](../references/neighborhood-validation.md) recipe for another representation or k, and calculate batch association across PCs and within-group mixing as appropriate. PCA axes show variance and permit loadings; nonlinear axes do not. Cluster distances, density, shapes and apparent trajectories in 2D do not establish relationships; validate against independent high-dimensional evidence. Record method, preprocessing, seed, parameters and limitations.

Done when the validation metric, representation and unresolved uncertainty are reported without treating the display as proof.
