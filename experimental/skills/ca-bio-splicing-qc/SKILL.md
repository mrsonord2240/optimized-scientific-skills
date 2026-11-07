---
name: ca-bio-splicing-qc
description: Use when assessing RNA-seq design and data quality for alternative splicing analysis.
category: Data Analysis
license: MIT
author: GPTomics
---

# Splicing Quality Control

Choose one route, read it, and run its command before writing new code. Paths are relative to this Skill.

| Need | Read |
|---|---|
| Junction saturation, annotation and support | `routes/junctions.md` |
| Cohort STAR two-pass alignment | `routes/star.md` |
| Strandedness, gene-body bias or rRNA | `routes/library.md` |
| Donor and acceptor sequence strength | `routes/sites.md` |
| Planned design or annotation choice | `routes/design.md` |
| Related analysis | `routes/related.md` |

Junction annotation and saturation require BED12 and matching BAM contigs. Thresholds are assay-dependent research QC, not patient variant classification.
Source checked: RSeQC 5.0.5, STAR 2.7.11b, pysam 0.24.1.
