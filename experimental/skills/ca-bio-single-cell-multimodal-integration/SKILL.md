---
name: ca-bio-single-cell-multimodal-integration
category: Data Analysis
description: Use when jointly analyzing multimodal single-cell data (CITE-seq RNA plus protein, 10x Multiome RNA plus ATAC, or unpaired and mosaic RNA and ATAC) and choosing the joint method.
tool_type: mixed
primary_tool: Seurat
license: MIT
author: GPTomics
---

# Multimodal Single-Cell Integration

Pick the one row that matches your data. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your data and task | Read |
|--------------------|------|
| Not sure which method fits the design (paired, unpaired, mosaic) | `routes/choose.md` |
| CITE-seq (RNA + ADT), Seurat: DSB denoising, then WNN | `routes/cite-seq-wnn-r.md` |
| CITE-seq, Python only, no denoising (CLR + Muon WNN) | `routes/cite-seq-wnn-python.md` |
| CITE-seq, denoised protein and foreground probabilities (totalVI) | `routes/totalvi.md` |
| Multiome (RNA + ATAC, same cell), Seurat/Signac WNN | `routes/multiome-wnn.md` |
| Shared and modality-specific factors (MOFA+) | `routes/mofa.md` |
| Mosaic: some cells RNA + ATAC, some one modality (MultiVI) | `routes/multivi.md` |
| Unpaired scRNA and scATAC, no shared cells (GLUE) | `routes/glue.md` |
| Unpaired scATAC onto an scRNA reference through a multiome bridge | `routes/bridge.md` |
| Wrong-looking result, MuData housekeeping, one modality dominates | `routes/errors.md` |
| Per-modality QC, clustering, batch correction, annotation, chromVAR, enrichment | `routes/other-skills.md` |

Rules on every route:

- Classify the task by anchor structure first (paired, unpaired, mosaic); it decides which methods apply.
- Denoise each modality in its native pipeline before the joint embedding. ADT background is ambient antibody, cell-intrinsic non-specific binding and spillover: CLR alone does not remove it. Use DSB (needs the raw matrix with empty droplets) or totalVI.
- Imputed modalities are inferences, not measurements. Flag any DE or marker call on MultiVI, StabMap, Cobolt or gene-activity values as model-dependent.
- Report the per-cell modality weights from WNN and check clustering survives down-weighting a dominant modality.
- A shared genome build and one common peak set across all ATAC cells and datasets.
- Seed model training (scvi-tools and GLUE are stochastic); the scripts do.
- RNA and protein disagreeing on a marker is often real biology; do not correct it away.

Tested with scanpy 1.10, Seurat 5.0 (bridge on 5.5.0 with Signac 1.17.1), anndata 0.10, scvi-tools 1.5.1, mudata 0.4.1, dsb 2.0.1.
