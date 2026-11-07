---
name: ca-bio-single-cell-clustering
category: Data Analysis
description: Use when clustering preprocessed single-cell RNA-seq cells into groups, choosing a resolution, or judging whether clusters are real.
tool_type: mixed
primary_tool: Seurat
license: MIT
author: GPTomics
---

# Single-Cell Clustering

Pick the one row that matches your data. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your data and task | Read |
|--------------------|------|
| Preprocessed AnnData (`.h5ad`); cluster, sweep resolution, UMAP | `routes/scanpy.md` |
| Preprocessed Seurat object (`.rds`); cluster, sweep resolution, UMAP | `routes/seurat.md` |
| A chosen partition to check for batch confounding and stability | `routes/validate.md` |
| Resolve fine states inside one cluster | `routes/subcluster.md` |
| A formal test that a split is real, or marker p-values offered as proof | `routes/significance.md` |
| Something looks wrong (one blob, too many clusters, labels change) | `routes/errors.md` |

Rules on every route:

- Input is normalized, log-transformed data with highly variable genes chosen. Raw counts are not clustering input.
- Resolution is a scale, not a truth. Sweep it, report cluster counts and sizes per value, and pick the coarsest level that orthogonal evidence supports. Do not tune it to match known cell types.
- `n_pcs` matters far more than `n_neighbors`. The elbow is a lower bound: check stability across nearby values, usually in 30-50.
- Cluster on the graph or PCA, never on UMAP or tSNE coordinates; embedding distances, sizes and gaps are not metric.
- Marker p-values after clustering rank genes only. They are not evidence that clusters differ (the same data chose the labels).
- Pin the Leiden backend and seed so labels reproduce across versions.

Tested with scanpy 1.12.4, anndata 0.13.3, igraph 1.0.0, leidenalg 0.12.0, Seurat 5.5.0, clustree.
