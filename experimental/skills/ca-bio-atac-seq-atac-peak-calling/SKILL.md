---
name: ca-bio-atac-seq-atac-peak-calling
description: Use when calling accessible chromatin peaks from aligned ATAC-seq reads.
category: Data Analysis
license: MIT
author: GPTomics
---

# ATAC-seq peak calling

Choose and read one route, then run its command before writing code. Paths are relative to this skill.

| Need | Route |
|---|---|
| Two replicates and reproducibility by IDR | `routes/replicates.md` |
| Single sample, fragment or broad peaks | `routes/single.md` |
| Joint calling or nucleosome-aware HMM | `routes/alternative.md` |

Use deduplicated, MAPQ-filtered, mitochondrial-read-removed BAMs and an assembly-matched blacklist. Select read-length-matched effective genome size; retain input and intermediate files. For upstream QC use `bio-atac-seq-atac-qc` and for consensus summits use `bio-atac-seq-consensus-peakset`. GPTomics bioSkills, MIT; notice in `LICENSE`; upstream d91ed3d563019e649dc854c56ccd62551359488a.

Source checked MACS3 3.0.4, samtools 1.24, bedtools 2.31.1 and IDR 2.0.4.2.
