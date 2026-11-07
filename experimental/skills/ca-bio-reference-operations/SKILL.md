---
name: ca-bio-reference-operations
category: Data Analysis
description: Use when preparing reference FASTA files, resolving contig names, or generating consensus from alignments.
license: MIT
author: GPTomics
---

# Reference operations

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Need | Route |
|---|---|
| Prepare FASTA index, dictionary, sizes | `routes/prepare.md` |
| Fetch or subset reference sequence | `routes/fetch.md` |
| Consensus from alignments, or compare to reference | `routes/consensus.md` |
| Resolve or rename contigs | `routes/contigs.md` |
| CRAM reference resolution and related work | `routes/related.md` |

Match reference sequence identity and contig names to the alignment before using either.

Tested with samtools 1.24, bcftools 1.24, pysam 0.24.1, GATK 4.6.2.0 and Picard 3.5.0.
