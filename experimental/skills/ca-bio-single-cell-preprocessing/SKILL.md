---
name: ca-bio-single-cell-preprocessing
category: Data Analysis
description: Use when quality-controlling, normalizing and selecting features for single-cell RNA-seq counts before dimensionality reduction and clustering.
tool_type: mixed
primary_tool: Seurat
license: MIT
author: GPTomics
---

# Single-Cell Preprocessing

Pick the one row that matches your data. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your data and task | Read |
|--------------------|------|
| 10x output or `.h5ad` in Python: QC, normalize, HVGs, PCA | `routes/scanpy.md` |
| 10x output in R / Seurat: QC, normalize, variable features | `routes/seurat.md` |
| Off-target markers everywhere (hemoglobin, immunoglobulin); ambient RNA | `routes/ambient.md` |
| Choosing a normalization, an HVG flavor, or whether to scale or regress out | `routes/normalization.md` |
| Something looks wrong (cell type vanished, garbage HVGs, almost nothing survives QC) | `routes/errors.md` |

Rules on every route:

- Order: raw matrix, empty-droplet calling, ambient removal (optional), cell and gene QC, doublet detection, normalization, HVGs, optional scaling, PCA. Ambient correction needs the raw matrix and comes before QC.
- Steps up to doublet detection are per sample, before any merge. Merging first leaks batch effects into every threshold.
- QC thresholds are MAD-adaptive (5 MAD on log counts, log genes and top-20-genes fraction; 3 MAD on mito fraction), with a tissue-aware mito cap. Never copy a PBMC cap into other tissue; for nuclei or unknown tissue use MAD only.
- Normalize raw counts exactly once. Keep integer counts in `layers['counts']`.
- Feed each HVG flavor its input: `seurat_v3` and Pearson residuals want raw counts; `seurat` and `cell_ranger` want log-normalized data. A wrong input runs silently and gives garbage.
- Do not reflexively scale or regress out mito or total counts; both are confounded with real cell state.
- Mito `^MT-` matches gene symbols only. With Ensembl IDs it matches nothing and the mito filter does nothing.

Tested with scanpy 1.12.4, scikit-misc 0.5.2, Seurat 5.5.0, SoupX 1.6.2, scran 1.34.0, celda 1.22.1.
