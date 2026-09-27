# Variant Normalization Usage Guide

## Overview

Normalize VCF files to a canonical representation (left-aligned, parsimonious, split, and optionally atomized) before comparing callers, merging cohorts, or matching against annotation databases.

## Prerequisites

- bcftools installed (`conda install -c bioconda bcftools`)
- Reference FASTA file (same one used for variant calling), indexed (`samtools faidx reference.fa`)
- Python parsing: `cyvcf2` (`pip install cyvcf2`)

## Quick Start

Tell your AI agent what you want to do:
- "Normalize my VCF by left-aligning indels and splitting multiallelic sites"
- "Prepare my VCF for annotation by normalizing variants against the reference"
- "Compare variants from two different callers after normalizing both VCFs"
- "Split multiallelic sites but keep SNPs and indels separate"

## Example Prompts

### Representation and Cross-Tool Decisions

> "My two cohorts were normalized with different tools and now show extra private variants -- reconcile their variant representation"

> "Decompose MNPs to match dbSNP but keep a haplotype-resolved copy for functional annotation"

> "A pathogenic ClinVar variant is showing as absent after annotation -- check whether my indels are left-aligned against the right reference"

> "The HGVS c. position my pipeline reports does not match the VCF POS for this deletion -- explain whether that is correct"

> "Split a multiallelic VCF while keeping AD and PL correct per allele"

## What the Agent Will Do

1. Check input VCF sorting and index status
2. Run bcftools norm with left-alignment and multiallelic splitting as needed
3. Remove duplicate records introduced by splitting
4. Verify normalization results by comparing variant counts before and after
5. Index the output VCF for downstream use

The pipeline order, tool commands, split/join/atomize options, and the HGVS clash are in
`SKILL.md`. `examples/normalize_vcf.sh` is a ready-made pipeline script; `examples/check_normalization.py`
estimates how many variants in a VCF need normalization.

## Tips

- Always normalize before comparing VCFs from different callers -- indel representation varies
- Standardize on ONE normalization tool and exact flag set across every cohort intended for comparison; `vt decompose_blocksub` splits MNPs by default while `bcftools norm` needs `--atomize`, so mixing tools manufactures spurious private variants
- Normalize against the EXACT reference used downstream (the annotation-database build) -- a one-base-off indel in a homopolymer silently misses its ClinVar/dbSNP entry with no error thrown
- Decompose for variant matching, but annotate functional consequence on the un-atomized (haplotype-resolved) VCF -- atomized adjacent SNVs in one codon can be misannotated
- Never hand-derive HGVS from POS: VCF left-aligns 5' on the genome, HGVS uses the transcript 3'-rule, so the two positions legitimately differ
- Normalization can increase variant count because multiallelic sites become multiple records

## Related Skills

- variant-calling/vcf-manipulation - Compare normalized VCFs with bcftools isec
- variant-calling/variant-annotation - Annotate after normalization
- variant-calling/filtering-best-practices - Filter after normalization
- variant-calling/vcf-basics - View and query normalized VCF files
