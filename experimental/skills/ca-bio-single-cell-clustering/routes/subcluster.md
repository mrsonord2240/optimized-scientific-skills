# Subclustering one cluster

```bash
python scripts/subcluster.py clustered.h5ad sub.h5ad cluster_key=leiden cluster=3 resolution=0.4
```

Subsets the cluster, then recomputes highly variable genes (2000), scaling, PCA (30), the kNN graph (20 PCs) and Leiden on the subset only, with the backend pinned. `X` must be log-normalized, not scaled. Writes `sub.h5ad` with a `leiden_sub` column.

- Never reuse the global PCA: its axes are uninformative inside a homogeneous subset and manufacture artificial sub-splits.
- Subclustering compounds double-dipping (the cells were selected by a clustering already). Do not use distinct markers or their p-values as a stop rule.
- Stop, and merge or keep the split only as a hypothesis, when it aligns with batch, sample or QC, fails resampling stability (`routes/validate.md`), fails a formal split test (`routes/significance.md`), or does not reproduce in an independent sample.
- A resolution can always split further; that is not a reason to continue.

Done when `sub.h5ad` exists and each sub-split has a keep or merge call with its evidence.
