# Psix trajectory

```bash
python scripts/psix_trajectory.py cells.h5ad psi_matrix.tsv mrna_matrix.tsv psix_results.tsv
```

Inputs are observed PSI and mRNA event-by-cell matrices and an h5ad containing X_pca for matching cells. The wrapper uses up to the first 20 PCA dimensions and four jobs, preserving the source example; `--jobs` changes parallelism. Interpret qvals and document any changes to method defaults; see `references/psix.md`. Do not impute PSI.

Done when trajectory-associated events and support are reported.
