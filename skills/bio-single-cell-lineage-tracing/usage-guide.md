# Lineage Tracing - Usage Guide

## Overview

Lineage tracing reads heritable marks (CRISPR/Cas9 scars, static expressed barcodes, or somatic mtDNA mutations) across single cells and asks which cells share which marks to recover an ontogenetic phylogeny or clonal grouping. The central idea: transcriptomic state does not fully predict fate (Weinreb 2020), so lineage information is orthogonal to expression and reconstruction is phylogenetics on error-prone scars where homoplasy and dropout dominate the result.

## Prerequisites

Cassiopeia and CoSpar/scanpy need two separate environments (numpy version conflict) —
see SKILL.md's "Installation and Version Compatibility" section for the exact commands.

## Quick Start

Tell your AI agent what you want to do:
- "Build a lineage tree from my CRISPR scar character matrix"
- "Compare solvers and tell me how robust the topology is"
- "Group clonally related cells from mtDNA mutations"
- "Recover early fate bias from sparse clonal barcodes"

## Example Prompts

### Tree Reconstruction
> "Reconstruct a lineage tree with Cassiopeia and collapse mutationless edges"
> "My data has heavy homoplasy and dropout - which solver should I use?"
> "Run a panel of solvers and compare them with Robinson-Foulds and triplets-correct"

### Clonal Dynamics and State
> "Integrate clones with transcriptomic state using CoSpar"
> "Compute early fate bias for monocyte vs neutrophil from clonal data"
> "Which clones expanded over time and what is their signature?"

### Mitochondrial and Native Tracing
> "Use mtDNA mutations to group clonally related cells in this human sample"
> "Why is my mtDNA giving clonal blobs instead of a deep tree?"

### Lineage vs State
> "Does transcriptomic state predict fate in my system, or is lineage adding new information?"
> "Validate my state-based trajectory against the lineage tree"

## What the Agent Will Do

1. Choose the recording assay by the question (deep topology, clonal state->fate, or retrospective native tissue)
2. Build a character matrix keeping missing (-1) distinct from unedited (0)
3. Assess missingness, informativeness, and barcode-collision risk before reconstruction
4. Reconstruct with a solver matched to scale and to homoplasy/dropout severity
5. Weight indels by formation probability and run a panel of solvers for robustness
6. Compare trees with Robinson-Foulds and depth-stratified triplets-correct
7. Integrate clone with state (CoSpar) to recover hidden fate bias and test state->fate

See SKILL.md's Governing Principle, Assay Decision Table, and Common Errors sections for
the underlying reasoning (state vs. fate, missing-vs-unedited, homoplasy, solver choice,
mtDNA clonal grouping) — those facts live there, not here.

## Related Skills

single-cell/trajectory-inference - state-based pseudotime/velocity that lineage data tests and corrects
single-cell/preprocessing - QC, doublet handling, and normalization upstream of barcode and clone calls
single-cell/clustering - cell-type labels annotated onto tree leaves and clones
phylogenetics/modern-tree-inference - general phylogenetic inference, parsimony vs ML, and branch support
