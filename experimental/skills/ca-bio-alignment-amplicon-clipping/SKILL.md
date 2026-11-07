---
name: ca-bio-alignment-amplicon-clipping
category: Data Analysis
description: Use when clipping PCR primer footprints from aligned amplicon reads.
license: MIT
author: GPTomics
---

# Amplicon primer clipping

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Task | Route |
|---|---|
| Clip and repair an amplicon BAM | `routes/clip.md` |
| Use iVar or hard clipping | `routes/alternatives.md` |
| Duplicate handling and downstream work | `routes/related.md` |

Use the primer scheme and reference build used for alignment. Standard coordinate duplicate marking is inappropriate for amplicon reads.

Tested with samtools 1.24, pysam 0.24.1, iVar 1.4.4.
