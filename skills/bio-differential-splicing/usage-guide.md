# Differential Splicing - Usage Guide

## Overview
Detect alternative splicing changes between conditions with rMATS-turbo, leafcutter, MAJIQ V3, SUPPA2 or Shiba. Which tool is right depends on the design (replicate count, cohort heterogeneity, annotation availability); `SKILL.md` has the decision tree, install notes, per-tool commands and failure modes.

## Quick Start
Tell your AI agent what you want to do:
- "Find differential splicing between tumor and normal samples"
- "Run rMATS and leafcutter in parallel and report concordant hits"
- "Use MAJIQ-HET for differential splicing across a heterogeneous patient cohort"
- "Compare splicing between treatment and control with low replicate count"
- "Prioritize differential splicing events by combined statistical and biological significance"

## Example Prompts

### Standard Replicate Designs
> "I have n=3 vs n=3 RNA-seq BAMs; run rMATS-turbo with FDR<0.05 and |dPSI|>0.10, then require inclusion plus skipping coverage >=10 in at least half the replicates of each group."

> "Use leafcutter Dirichlet-multinomial GLM on intron clusters from regtools junctions for annotation-free differential splicing; include batch as a groups-file covariate if it is not aliased with condition."

### Heterogeneous Cohorts
> "I have 30 tumor and 30 normal samples from heterogeneous patients; use MAJIQ V3 HET module with posterior threshold P(|dPSI|>0.2)>0.95."

### Low Replicate Count
> "n=2 vs n=2 design - use leafcutter or Shiba."

### Result Prioritization
> "Compute combined score (-log10(FDR) * |dPSI|) and pull the top 50 events with NMD-status and protein-domain annotation."

> "Cross-reference top hits with eCLIP RBP target databases to identify candidate trans-regulators."

## Related Skills

See the Related Skills section of `SKILL.md`.
