# MSA Parsing - Usage Guide

## Overview

This skill focuses on parsing and analyzing multiple sequence alignments (MSAs). It covers extracting information, analyzing gaps and conservation, filtering sequences, and preparing alignments for downstream analysis like phylogenetics or structure prediction. Install notes, thresholds, pitfalls and all code live in `SKILL.md`.

## Quick Start

Tell your AI agent what you want to do:
- "Find all the fully conserved positions in this alignment"
- "Remove sequences with more than 10% gaps"
- "Generate a consensus sequence from this alignment"

## Example Prompts

### Analyzing Content
> "Show me the composition of each column in the alignment"

> "Find positions that are conserved in at least 80% of sequences"

> "Count the gaps in each sequence"

### Filtering and Cleaning
> "Remove columns with more than 50% gaps"

> "Filter out sequences with too many gaps"

> "Remove duplicate sequences from the alignment"

### Extracting Information
> "Get the sequence for species_A from this alignment"

> "Extract columns 100-200 from the alignment"

> "List all sequence IDs in the alignment"

### Consensus and Conservation
> "Generate a consensus sequence with 70% threshold"

> "Find the most conserved regions in this alignment"

> "What is the consensus at each position?"

### Weighting, Coevolution and Reliability
> "Weight the sequences and report Neff for this alignment"

> "Which columns of this alignment are unreliable? Mask them"

> "Look for coevolving column pairs, and tell me whether the alignment is deep enough"

## Related Skills

See the Related Skills section of `SKILL.md`.
