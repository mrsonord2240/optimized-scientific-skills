# Alignment Statistics - Usage Guide

## Overview
Generate QC statistics from alignment files including mapping rates, read counts, coverage depth, and per-chromosome distributions. Commands, recipes, thresholds and caveats live in `SKILL.md`; install notes are under "Version Compatibility" there.

## Quick Start
Tell your AI agent what you want to do:
- "Get quick statistics for my BAM file"
- "Calculate coverage depth across my genome"
- "Generate a comprehensive QC report with plots"
- "Check per-chromosome read distribution"

## Example Prompts

### Quick Statistics
> "Run flagstat on my BAM file and interpret the results"

> "How many reads are mapped in sample.bam?"

> "What percentage of reads are properly paired?"

### Per-Chromosome Analysis
> "Show read counts per chromosome using idxstats"

> "Calculate the mitochondrial contamination percentage"

> "Check X/Y ratio for sex determination"

### Coverage Analysis
> "Calculate mean coverage depth for my WGS sample"

> "What percentage of the genome is covered at 20x?"

> "Generate depth statistics for target regions"

### Comprehensive Reports
> "Generate samtools stats and create QC plots"

> "Create a summary table of statistics for all my samples"

## Related Skills
See the "Related Skills" list at the end of `SKILL.md`.
