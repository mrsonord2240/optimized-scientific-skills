# Differential Abundance Testing - Usage Guide

## Overview

Differential abundance testing asks whether cell-type proportions or composition changed between conditions (treatment vs control, disease vs healthy). Because proportions live on a simplex (they sum to 1 and are not independent), naive per-cluster proportion tests are invalid. This skill covers cluster-free neighborhood testing (Milo) and cluster-based compositional models (scCODA, sccomp, propeller), and how to keep compositional shifts from being misread as differential expression. Install commands, method choice, thresholds and failure modes are in `SKILL.md`.

## Quick Start

Tell your AI agent what you want to do:
- "Did any cell-type proportions change between my conditions?"
- "Test differential abundance with Milo on my integrated data"
- "Run scCODA to find populations that expanded with treatment"
- "Check whether my DE signal is actually a compositional shift"

## Example Prompts

### Cluster-free abundance
> "Run Milo differential abundance on my kNN graph and report SpatialFDR neighborhoods"
> "Find transitional states that expanded with treatment without committing to clusters"

### Cluster-based composition
> "Test cell-type proportion changes with scCODA using a stable reference cell type"
> "Run sccomp on my count table, robust to outlier samples"
> "Use propeller to quickly test proportion differences across groups"

### Guarding the DE/DA confound
> "Pair my pseudobulk DE with a differential-abundance test"
> "Is the change in this cluster's expression real or just a shift in substate proportions?"

### Interpretation
> "Which population actually drives the compositional change relative to the reference?"
> "Do Milo and scCODA agree on which cell types changed?"

## Related Skills

- clustering - Define the clusters whose abundance is tested
- cell-annotation - Annotate cell types before testing their proportions
- markers-annotation - Pair condition DE with abundance testing to separate the confound
- batch-integration - Build the integrated embedding Milo's kNN graph relies on
- differential-expression/deseq2-basics - Pseudobulk condition DE that abundance testing complements
- pathway-analysis/go-enrichment - Characterize the populations that expanded or contracted
