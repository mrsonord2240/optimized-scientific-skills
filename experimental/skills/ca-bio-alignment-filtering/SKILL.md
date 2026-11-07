---
name: ca-bio-alignment-filtering
category: Data Analysis
description: Use when filtering aligned reads by flags, mapping quality, regions, tags, or subsampling.
license: MIT
author: GPTomics
---

# Alignment filtering

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Need | Route |
|---|---|
| Flags, MAPQ, assay-aware selection | `routes/filter.md` |
| Region or BED extraction | `routes/region.md` |
| Tags, expressions, read groups | `routes/expression.md` |
| Fraction or target read count | `routes/subsample.md` |
| Python filtering | `routes/python.md` |

Read-level filters can orphan mates. Check counts before and after and validate output. Duplicate removal requires duplicates already marked; avoid coordinate marking on amplicons.

Tested with samtools 1.24 and pysam 0.24.1.
