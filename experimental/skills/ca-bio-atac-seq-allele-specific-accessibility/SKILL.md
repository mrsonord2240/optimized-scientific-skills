---
name: ca-bio-atac-seq-allele-specific-accessibility
category: Data Analysis
description: Use when analyzing allele-specific accessibility from ATAC-seq and heterozygous genotypes.
license: MIT
author: GPTomics
---

# Allele-specific accessibility

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Need | Route |
|---|---|
| Individual phased peak analysis | `routes/individual.md` |
| Moderate cohort joint analysis | `routes/rasqual.md` |
| Other designs or no genotypes | `routes/other.md` |

Match sample, assembly, contigs and reference; apply WASP before counting and verify the null reference-allele distribution. Pool only phased alleles within a phase set. Report depth, effect, multiplicity, mapping-bias correction and limitations.

Tested contract: WASP d3b8447, GATK 4.6, samtools/bcftools 1.21, bowtie2 2.5.
