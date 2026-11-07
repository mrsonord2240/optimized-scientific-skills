---
name: ca-bio-long-read-splicing
description: Use when discovering, quantifying, or testing alternative splicing and isoforms from long-read RNA sequencing.
category: Data Analysis
license: MIT
author: GPTomics
---

# Long-Read Splicing

Choose a route, read it, and run its command before writing new code. Paths are relative to this Skill.

| Data or question | Read |
|---|---|
| End-to-end annotated HiFi or ONT analysis | `routes/pipeline.md` |
| Platform-aware alignment or microexons | `routes/alignment.md` |
| FLAIR correction, collapse and differential splicing | `routes/flair.md` |
| IsoQuant or Bambu discovery and quantification | `routes/discovery.md` |
| SQANTI3 classification and filtering | `routes/sqanti.md` |
| rMATS-long differential events | `routes/rmats.md` |
| Differential transcript usage from long-read counts | `routes/dtu.md` |
| Kinnex single-cell library processing | `routes/single-cell.md` |
| Related analyses | `routes/related.md` |

Review intra-priming and RT-switching before reporting novel isoforms. Use Bambu NDR 0.1 initially and stop novel claims on an NDR-approximated warning.
Source checked versions: minimap2 2.31, FLAIR 3.0.1, IsoQuant 4.0.0, Bambu 3.8.3, SQANTI3 6.0.2, rMATS-long 2.1.0.
