---
name: ca-bio-pileup-generation
description: Use when generating alignment pileups, counting alleles at genomic positions, or calling variants from pileups.
category: Data Analysis
license: MIT
author: GPTomics
---

# Pileup Generation

Choose one route, read it, and run its command before writing code. Paths are relative to this Skill.

| Need | Read |
|---|---|
| Text pileup for positions, regions or samples | `routes/text.md` |
| Allele counts or custom position analysis in Python | `routes/alleles.md` |
| Pileup-based variant calls or BCF intermediate | `routes/calling.md` |
| Related analysis | `routes/related.md` |

Match reference contigs to BAM headers; a missing reference contig can produce N even with exit status zero.
Checked source versions: samtools 1.24, bcftools 1.24, pysam 0.24.1.
