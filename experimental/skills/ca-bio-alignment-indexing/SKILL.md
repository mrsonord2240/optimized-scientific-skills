---
name: ca-bio-alignment-indexing
category: Data Analysis
description: Use when indexing alignment files or retrieving indexed genomic regions and statistics.
license: MIT
author: GPTomics
---

# Alignment indexing

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Need | Route |
|---|---|
| Create or refresh BAM or CRAM indices | `routes/index.md` |
| Fetch regions or counts | `routes/fetch.md` |
| Per-contig index statistics | `routes/stats.md` |
| Index or retrieve FASTA | `routes/fasta.md` |

Indexing requires coordinate-sorted alignments. Contig names must match exactly.

Tested with samtools 1.24 and pysam 0.24.1.
