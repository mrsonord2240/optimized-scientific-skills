---
name: ca-bio-duplicate-handling
category: Data Analysis
description: Use when marking or removing PCR and optical duplicates in an aligned DNA BAM, or deciding whether duplicate marking is wrong for the assay.
tool_type: cli
primary_tool: samtools
license: MIT
author: GPTomics
---

# Duplicate Handling

Pick the one row that matches your BAM. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your BAM | Read |
|----------|------|
| WGS, exome, capture, somatic, ChIP-seq, CUT&RUN, ATAC, ancient DNA; no UMIs | `routes/markdup.md` |
| Same, but you want it in Python | `routes/pysam.md` |
| UMIs (RX tag, or CB+UB tags from 10x), ctDNA, duplex | `routes/umi.md` |
| Bulk RNA-seq, scRNA, amplicon or hotspot panel, native long reads, 16S/ITS | `routes/not-markdup.md` |
| samblaster, Picard, biobambam2, sambamba, PacBio HiFi `pbmarkdup`, mapDamage | `routes/other-markers.md` |
| A duplicate rate, or a BAM with the duplicates filtered out | `routes/rate-and-filter.md` |
| An error message | `routes/errors.md` |

Rules on every route:

- Mark, do not remove. `samtools markdup -r` is irreversible; filter on flag 1024 downstream instead.
- ChIP-seq and CUT&RUN: mark only, never remove; the peak caller handles duplicates (`macs3 --keep-dup auto`).
- ATAC-seq: mark before the Tn5 +4/-5 shift, not after.
- Do not hand-roll a marker in Python on (chrom, pos, strand): it ignores soft-clip positions, read groups and supplementary alignments.
- A 50-95% duplicate rate on a 10x or RNA BAM means the wrong tool, not a bug.

Tested with samtools 1.24, pysam 0.24.1, umi_tools 1.1.6, fgbio 4.1.1, Picard 3.5.0, samblaster 0.1.26, sambamba 1.0.1, pbmarkdup 1.2.0, mapDamage 2.2.2.
