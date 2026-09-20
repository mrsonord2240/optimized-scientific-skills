# Pairwise Alignment - Usage Guide

## Overview

This skill performs pairwise sequence alignment to compare two DNA, RNA, or protein sequences. It uses Biopython's `PairwiseAligner` class which implements dynamic programming algorithms for finding optimal alignments.

## Quick Start

Tell your AI agent what you want to do:

- "Align these two DNA sequences and show me the best alignment"
- "Compare this protein sequence against a reference using BLOSUM62"
- "Find the best matching region between these two sequences"

## Example Prompts

### Global Alignment
> "Perform a global alignment between ACCGGTAACGTAG and ACCGTTAACGAAG"

> "Align the first two sequences in my FASTA file"

### Local Alignment
> "Find the best local alignment between these sequences to identify conserved regions"

> "Use Smith-Waterman to find matching regions in these proteins"

### Protein Alignment
> "Align these two protein sequences using BLOSUM62 scoring"

> "Compare my query protein against the reference with appropriate gap penalties"

### Scoring and Analysis
> "Calculate the alignment score between these sequences"

> "Show me all optimal alignments and their scores"

## Details

Modes, scoring, gap conventions, significance testing and when to leave pairwise DP behind are in `SKILL.md`.
