# Reference Operations - Usage Guide

## Overview
Work with reference genomes including indexing FASTA files, extracting sequences, creating sequence dictionaries, and generating consensus sequences from alignments.

Commands, options, contig-naming tables, error messages and Python code live in `SKILL.md`; the guide only lists what to ask for. Install notes are under "Version Compatibility" there.

## Quick Start
Tell your AI agent what you want to do:
- "Index my reference FASTA file"
- "Extract sequence for chr1:1000-2000"
- "Create a sequence dictionary for GATK"
- "Generate consensus sequence from my alignments"

## Example Prompts

### FASTA Indexing
> "Index my reference genome for random access"

> "Prepare the reference for variant calling with GATK"

> "Create all required index files for my reference"

### Sequence Extraction
> "Extract the sequence for chr1:1000000-2000000"

> "Get multiple genomic regions and save to FASTA"

> "Extract the reverse complement of a region"

### Sequence Dictionary
> "Create a sequence dictionary for my reference"

> "Prepare reference files for Picard tools"

### Contig Names
> "My BAM says chr22 but my reference says 22 -- fix it without re-aligning"

> "Which GRCh38 reference is this BAM aligned to?"

### Consensus Generation
> "Generate a consensus sequence from my BAM file"

> "Create consensus for a specific region with minimum depth"

> "Compare consensus to reference to find differences"

## What the Agent Will Do

1. Check that the reference exists and is readable
2. Create the index files (FAI for samtools, `<name>.dict` for GATK/Picard)
3. Extract sequences or generate the consensus as requested
4. Verify the output files were created and match the reference
