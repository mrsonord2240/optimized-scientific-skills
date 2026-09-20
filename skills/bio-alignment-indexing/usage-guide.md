# Alignment Indexing - Usage Guide

## Overview
Create and use indices for random access to BAM and CRAM files, enabling fast region queries without reading entire files. Commands, index types, CRAM reference handling, staleness checks and error messages are in `SKILL.md` (sections: Index Types, CRAM, Index Staleness, Common Errors); install notes are under Version Compatibility there.

## Quick Start
Tell your AI agent what you want to do:
- "Index my BAM file for random access"
- "Get per-chromosome read counts from the index"
- "Check if my BAM file has an index, and whether it is stale"
- "Create a CSI index for large chromosomes"

## Example Prompts

### Creating Indices
> "Create an index for sample.bam"

> "Index every BAM in this folder, skipping ones that already have a fresh index"

> "Create CSI index for my genome with large chromosomes"

> "Index my CRAM and pull the reads for chr1:1000-2000 (reference is ref.fa)"

### Using Indices
> "Count reads per chromosome using idxstats"

> "Extract reads from chr1:1000000-2000000"

> "Get reads overlapping regions in my BED file"

### FASTA Indexing
> "Index my reference FASTA for random access"

> "Extract sequence for chr1:1000-2000 from reference"

### Index Statistics
> "Show per-chromosome read distribution"

> "Calculate mitochondrial contamination percentage"

## What the Agent Will Do

1. Verify the BAM file is coordinate-sorted (required for indexing)
2. Create the appropriate index type (BAI, CSI, or CRAI)
3. Place the index file alongside the BAM file
4. Verify the index was created successfully
5. Use the index for efficient region queries or statistics

## Related Skills
See `SKILL.md` (Related Skills).
