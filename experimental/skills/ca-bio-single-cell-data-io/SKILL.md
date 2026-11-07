---
name: ca-bio-single-cell-data-io
category: Data Analysis
description: Use when loading 10x Cell Ranger output or moving single-cell objects between AnnData, Seurat and SingleCellExperiment without losing data.
tool_type: mixed
primary_tool: Seurat
license: MIT
author: GPTomics
---

# Single-Cell Data I/O

Pick the one row that matches your task. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your task | Read |
|-----------|------|
| Load 10x Cell Ranger output (matrix directory or `.h5`), raw or filtered, into Python or R | `routes/load-10x.md` |
| Where counts, normalized values, metadata and embeddings go in AnnData or Seurat v5 | `routes/object-structure.md` |
| Convert between h5ad, Seurat and SingleCellExperiment; recover `.raw` | `routes/convert.md` |
| Pick a storage format (h5ad, zarr, RDS, h5mu, loom), or a read/write default surprised you | `routes/formats-and-defaults.md` |
| Something went missing or swapped (layers, cells and genes, ADT features, wrong values) | `routes/errors.md` |

Rules on every route:

- The dominant failure is silent: the matrix changes meaning during a read or conversion. Diff slot inventories before and after any cross-ecosystem hop (`scripts/h5ad_inventory.py`) and keep the original file.
- AnnData is cells x genes; Seurat and SingleCellExperiment are genes x cells. A faithful conversion transposes the matrix and swaps which axis the metadata annotates.
- Keep matrices sparse. Dense `.toarray()` on 500k cells x 30k genes is about 60 GB.
- Keep the RAW (unfiltered) Cell Ranger matrix. EmptyDrops, SoupX, CellBender and DecontX need it. The filtered matrix is cell-called, not ambient-corrected.
- Prefer stable Ensembl IDs (`gene_ids`) over symbols when joining datasets.
- Verify a converter is maintained and matches the installed versions before relying on it. SeuratDisk is abandoned and broken on Seurat v5.

Tested with scanpy 1.12.4, anndata 0.13.3, Seurat 5.5.0, SeuratObject 5.4.0, zellkonverter 1.16.0, schard 1.1.0 (anndataR needs R 4.5 and is not installed).
