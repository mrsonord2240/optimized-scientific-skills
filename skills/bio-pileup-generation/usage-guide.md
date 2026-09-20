# Pileup Generation - Usage Guide

## Overview
Generate pileup data showing all reads covering each genomic position for variant calling and position-level analysis. Commands, flags, defaults, pysam parameters and troubleshooting are in `SKILL.md`; this guide is for choosing the Skill and phrasing requests.

## Quick Start
Tell your AI agent what you want to do:
- "Generate pileup for variant calling"
- "Count alleles at a specific position"
- "Call variants using bcftools"
- "Check read support at position chr1:1000000"

## Example Prompts

### Basic Pileup
> "Generate text pileup for my BAM file"

> "Create pileup for chr1:1000000-2000000 only"

> "Generate pileup with quality filtering (MAPQ 20, baseQ 20)"

### Variant Calling
> "Call variants using bcftools mpileup and bcftools call"

> "Generate a BCF file so I can re-run variant calling without repeating the pileup"

> "Call variants from multiple samples together"

### Position Analysis
> "Count allele frequencies at position chr1:1000000"

> "Find all positions with alternative alleles above 10%"

> "Check read support for known variants"

## What the Agent Will Do

1. Verify the BAM file is indexed and reference is available
2. Apply quality filters (mapping quality, base quality)
3. Generate pileup for the requested regions
4. Parse pileup output to extract allele counts
5. Call variants or report position-level statistics
6. Output results in requested format (text, BCF, VCF)

## Related Skills
See the Related Skills list in `SKILL.md`.
