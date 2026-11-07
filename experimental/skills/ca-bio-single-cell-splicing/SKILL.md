---
name: ca-bio-single-cell-splicing
category: Data Analysis
description: Use when analyzing single-cell splicing or alternative polyadenylation with chemistry-aware methods.
license: MIT
author: GPTomics
---

# Single-cell splicing

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Chemistry or need | Route |
|---|---|
| Assess library suitability | `routes/chemistry.md` |
| Full-length plate PSI and cell types | `routes/marvel.md` |
| Per-cell PSI and covariates | `routes/brie2.md` |
| Plate junction clusters | `routes/scquint.md` |
| Cell-state splicing discovery | `routes/spliz.md` |
| Trajectory association | `routes/psix.md` |
| 10X 3-prime alternative polyadenylation | `routes/sierra.md` |
| Replicate-aware differential splicing | `routes/pseudobulk.md` |
| Full-length isoforms | `routes/related.md` |

Check chemistry first: 10X 3-prime is unsuitable for transcriptome-wide per-cell splicing. Do not impute PSI or treat cells as independent biological replicates. Match genome build and contig names; filter doublets. Nuclear intron retention does not necessarily indicate mature isoform regulation.

Tested with MARVEL 2.0.5, BRIE2 2.3.0, scQuint 0.3.3, Psix 0.10.8 and Sierra 0.99.27; SpliZ was not run.
