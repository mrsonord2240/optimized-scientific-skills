# Alignment I/O - Usage Guide

## Overview

This skill handles reading, writing, and converting multiple sequence alignment (MSA) files. It uses Biopython's `AlignIO` module (plus `pyhmmer` for Pfam-scale streaming) which provides a consistent interface for many alignment formats used in phylogenetics and sequence analysis.

## Quick Start

Tell your AI agent what you want to do:
- "Read this Clustal alignment file and show me the sequences"
- "Convert my PHYLIP alignment to FASTA format"
- "Parse the Stockholm file and extract sequence IDs"

## Example Prompts

### Reading Alignments
> "Load the alignment from alignment.aln and tell me how many sequences it has"

> "Parse all alignments from this Stockholm file"

> "Read the PHYLIP file and show me the first 50 columns"

### Writing Alignments
> "Save this alignment as a FASTA file"

> "Write the alignment to Nexus format for MrBayes"

> "Export to PHYLIP-relaxed format to preserve long sequence names"

### Format Conversion
> "Convert my Clustal alignment to PHYLIP format"

> "Convert all .aln files in this directory to FASTA"

> "Change the alignment from Stockholm to Nexus"

### Working with Alignment Data
> "Extract sequences 1-10 from the alignment"

> "Trim the alignment to columns 50-150"

> "Get all sequence IDs from the alignment"

## What the Agent Will Do

1. Load alignment file using appropriate format parser
2. Access sequences, IDs, and alignment properties
3. Perform requested operations (slice, filter, convert)
4. Write output in specified format
5. Report alignment statistics (length, sequence count)

## Notes

Install steps, the downstream-format table (which format for IQ-TREE, MrBayes, HMMER, PAML), PHYLIP naming pitfalls and annotation loss on conversion are all in `SKILL.md` ("Version Compatibility", "Format Selection for Downstream Tools", "Format-Specific Notes").
