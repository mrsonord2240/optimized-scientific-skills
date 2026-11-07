# Check a chosen partition

```bash
python scripts/validate_partition.py clustered.h5ad leiden_r0.6 covariates=batch,sample,condition resolution=0.6 n_bootstrap=20
```

The object must already hold PCA (`X_pca`) and the cluster column. The script prints the adjusted Rand index between the partition and each covariate present in `.obs`, then re-clusters 80% subsamples `n_bootstrap` times and prints each cluster's mean best-match Jaccard, lowest first. It writes `<input>.stability.tsv`.

- Reject or flag a split when it aligns with a covariate (high ARI) or has Jaccard below about 0.6-0.7.
- Passing is necessary, not sufficient. A reproducible cluster can still be cell cycle, dissociation stress (FOS, JUN, HSPA1A), mitochondrial fraction, ambient RNA or batch. Check that the cluster survives regressing cycle and mito out and carries non-stress markers.
- A claim that a split is a distinct population also needs `routes/significance.md` and an independent sample.
- Resolution must match the one that made the partition, or the stability numbers describe a different partition.

Done when ARI per covariate and Jaccard per cluster are reported, with a keep/flag/merge call per cluster.
