# Cluster with Scanpy

```bash
python scripts/cluster_scanpy.py preprocessed.h5ad clustered.h5ad n_pcs=30 n_neighbors=15 resolutions=0.2,0.4,0.6,0.8,1.0
```

Runs PCA (50 components, if `X_pca` is absent), builds the kNN graph on `n_pcs`, runs Leiden at every resolution with the backend pinned (`flavor='igraph', n_iterations=2, directed=False`, seeded), computes UMAP, and writes `clustered.h5ad` with one `leiden_r<res>` column per resolution plus `<out>.sweep.tsv` (cluster count and size range per resolution). It refuses input that looks like raw counts.

- Report the sweep table: cluster count and size distribution per resolution, the chosen coarsest level, and the evidence for it from `routes/validate.md`.
- Each cluster's batch, sample and QC composition belongs in the report. A cluster that maps to one sample is a batch effect: integrate first (batch-integration Skill).
- Doublets form fake intermediate clusters: run doublet detection before clustering.
- Marker genes are annotation and ranking evidence only.
- `min_dist` for UMAP affects display, never the partition.
- Unpinned backends change labels between Scanpy versions. Switching `flavor` or `n_iterations` changes labels.

Done when `clustered.h5ad` and the sweep table exist and you have reported counts per resolution.
