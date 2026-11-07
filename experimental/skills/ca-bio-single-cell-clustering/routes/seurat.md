# Cluster with Seurat

```bash
Rscript scripts/cluster_seurat.R preprocessed.rds clustered.rds dims=30 resolutions=0.2,0.4,0.6,0.8,1.0 algorithm=4
```

Runs PCA (50 components), `FindNeighbors` on `dims`, `FindClusters` at every resolution, `RunUMAP`, writes `clustered.rds`, a clustree plot `clustree.pdf` (cell flow across the sweep) and `<out>.sweep.tsv` (cluster count per resolution). The columns are `RNA_snn_res.<r>`.

- `algorithm=4` is Leiden (needs the `leidenbase` package; the script stops with a message if it is missing). `algorithm=1` is Louvain, Seurat's own default, which can return internally disconnected communities.
- Seurat 5 with `algorithm=4` uses `leidenbase` unless `leiden_method` is set; check `?FindClusters` for the installed version.
- Sweep, report counts and sizes per resolution, pick the coarsest defensible level, and label finer splits as hypotheses (see `routes/scanpy.md` for the report items).
- Do not repeat the run to compare engines.

Done when `clustered.rds`, `clustree.pdf` and the sweep table exist.
