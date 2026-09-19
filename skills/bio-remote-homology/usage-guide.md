# Remote Homology Usage Guide

## Overview

Detect distant homologs using profile and structure-aware methods. Encodes the modern landscape: PSI-BLAST and jackhmmer for iterative profile search; HHblits/HHsearch for profile-profile (PDB70, Pfam); MMseqs2 and DIAMOND as the modern blastp replacements (100-10,000x faster); Foldseek (van Kempen 2024) for structure-aware homology via the 3Di alphabet, with ProstT5 for sequence-only access to structural search.

Setup, the decision matrix, tool-by-tool details, code patterns, failure modes and common errors all live in `SKILL.md` -- this guide only covers when to reach for this skill and what to ask for.

## Quick Start

- "Find structural homologs of this protein via Foldseek against AlphaFoldDB Swiss-Prot subset"
- "Run jackhmmer 3 iterations against UniRef90; save the final HMM for downstream hmmsearch"
- "Search a protein query with MMseqs2 at -s 7.5 (HMMER-equivalent sensitivity) instead of BLAST"
- "Annotate domains with hmmscan against Pfam-A using --cut_ga calibrated thresholds"
- "DIAMOND --ultra-sensitive against UniRef90 for a metagenomic protein set"

## Example Prompts

### Foldseek for divergent structural homologs

> "I have a protein where blastp finds nothing significant. Predict the structure (or use ProstT5 for sequence-only) and search AlphaFoldDB with Foldseek easy-search. Report TM-scores and probability cutoff -- hits with prob > 0.9 are structurally confident."

### PSI-BLAST with drift protection

> "Run psiblast for 3 iterations against UniRef90 with -inclusion_ethresh 0.002 (stricter than the 0.005 default). Save the PSSM with -out_pssm and the included sequence set so I can audit for paralog contamination. Don't iterate to convergence -- 4+ iterations drift."

### MMseqs2 as PSI-BLAST replacement

> "Same iterative profile search as PSI-BLAST but use MMseqs2 with --num-iterations 3 and -s 7.5 (HMMER-equivalent sensitivity). It's 100x faster."

### Pfam domain annotation

> "Annotate domains of these 5,000 protein sequences using hmmscan against Pfam-A with --cut_ga. Use the calibrated gathering thresholds, not arbitrary E-value cutoffs."

### HHsearch against PDB70

> "Build an HHblits profile of this protein against UniRef30 (3 iterations), then hhsearch against PDB70 for the deepest possible structural homology to known PDB entries."

### DIAMOND for metagenomic scale

> "I have 1 million predicted ORFs from a metagenome. Use DIAMOND blastp --ultra-sensitive against UniRef90 with -p 32. blastp would take days; DIAMOND will finish in an hour."

## Related Skills

- blast-searches - Remote BLAST baseline
- local-blast - Local BLAST+ for moderate-scale work
- ortholog-inference - Distinct topic: orthology vs homology (RBH, OrthoFinder, OMA)
- alignment/multiple-alignment - Build MSAs for HMM profiles
- structural-biology/alphafold-predictions - Predict structures for Foldseek queries
