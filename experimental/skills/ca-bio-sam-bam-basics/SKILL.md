---
name: ca-bio-sam-bam-basics
category: Data Analysis
description: Use when inspecting or converting SAM, BAM, and CRAM alignment files.
license: MIT
author: GPTomics
---

# SAM, BAM and CRAM basics

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Need | Route |
|---|---|
| View, count or inspect alignments | `routes/view.md` |
| Convert formats | `routes/convert.md` |
| Fetch one or multiple regions | `routes/regions.md` |
| Interpret flags, CIGAR, coordinates, MAPQ or tags | `routes/interpret.md` |
| Python or R access | `routes/programmatic.md` |

Use the exact reference for CRAM and confirm contig names. Preserve supplementary alignments for split-read analyses.

Tested with samtools 1.24, pysam 0.24.1 and bcftools 1.24.
