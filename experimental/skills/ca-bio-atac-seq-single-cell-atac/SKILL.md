---
name: ca-bio-atac-seq-single-cell-atac
description: Use when processing, clustering, or integrating single-cell ATAC-seq data.
category: Data Analysis
license: MIT
author: GPTomics
---

# Single-cell ATAC-seq

Choose and read one route, then run its command before writing code. Paths are relative to this skill.

| Need | Route |
|---|---|
| 10X ATAC or Multiome peaks in Signac | `routes/signac.md` |
| ArchR or SnapATAC2 workflow | `routes/ecosystems.md` |
| Doublets, cluster peaks or pseudobulk | `routes/peaks.md` |
| Paired RNA/ATAC integration or specialist analysis | `routes/multiome.md` |

Confirm chemistry, barcode identity, assembly, fragments index and sample metadata. Gene activity approximates expression; it is not direct RNA measurement. GPTomics bioSkills, MIT notice in `LICENSE`, upstream d91ed3d563019e649dc854c56ccd62551359488a.

Source checked Signac 1.17.1, Seurat 5.5.1, ArchR 1.0.3, SnapATAC2 2.10.0 and MACS3 3.0.4.
