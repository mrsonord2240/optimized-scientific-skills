# Alignment Amplicon Clipping

Trim PCR primers from aligned reads in amplicon-panel BAMs.

## Overview

Amplicon panels (SARS-CoV-2 ARTIC, hereditary cancer panels, ctDNA hot-spot panels, fusion panels, 16S rRNA) use designed PCR primers for enrichment. The primer footprint at read ends does not represent biological sequence; without trimming, downstream variant calls suppress true variants under primers and falsely confirm reference at primer locations.

This skill covers post-alignment primer clipping with `samtools ampliconclip`, choosing the clip mode for your read geometry, tag repair, and the relationship to duplicate marking. Commands, decision points, thresholds and failure modes are all in `SKILL.md`.

You need: samtools 1.11 or later, a coordinate-sorted indexed amplicon BAM, and the primer BED for your kit with strand in column 6 (install notes in `SKILL.md`).

## Example Prompts

### SARS-CoV-2 / ARTIC

> "I have a SARS-CoV-2 ARTIC BAM aligned with bwa-mem. Soft-clip primers using artic_v3 primer BED, then prepare for consensus calling."

> "Compare samtools ampliconclip vs iVar trim for the ARTIC SARS-CoV-2 pipeline."

> "My nanopore ARTIC BAM still has primer sequence at the 3' end of reads after ampliconclip."

### Cancer Panels

> "ctDNA hot-spot panel with TruSeq UMIs. Trim primers, then run UMI consensus rather than markdup."

> "Hereditary cancer panel from a Roche AVENIO design. Why is my variant caller missing variants at primer-overlapping positions?"

### General Amplicon Workflows

> "My amplicon BAM ran through markdup and now flagstat shows 99% duplicates. What went wrong?"

> "After ampliconclip, IGV shows no mismatches and my NM filter drops reads. What did clipping remove?"

> "Hard-clip primers for archival storage."

## Related Skills

Listed in `SKILL.md`.
