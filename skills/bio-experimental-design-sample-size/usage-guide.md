# Sample Size Estimation Usage Guide

## Overview

This guide covers estimating the minimum number of biological replicates (or cells and donors) for a target power at a target false discovery rate across genomics assays: bulk RNA-seq, scRNA-seq, ChIP-seq, ATAC-seq, methylation, and proteomics. The central principle is that sample size counts biological units, not measurements: technical replicates reduce measurement noise but add no degrees of freedom for biological inference, and for single-cell studies the number of donors, not cells, sets the power for population-level differential expression. The guide stresses estimating dispersion from pilot data rather than guessing a coefficient of variation, and treats the familiar "n=3" as a publication convention rather than a calculation. See `SKILL.md` for setup, code, thresholds, and failure modes.

## Quick Start

Tell your AI agent what you want to do:
- "How many biological replicates for bulk RNA-seq to detect 1.5-fold at 80% power, FDR 0.05?"
- "Estimate dispersion from my pilot and size the full study"
- "Is my budget-fixed n of 6 adequate, and what FDR does it actually achieve?"
- "How many donors versus cells for a scRNA-seq disease-versus-control study?"
- "Does my n=3 plan have enough power, or do I need more?"

## Example Prompts

### Bulk RNA-seq

> "I'm comparing tumor versus normal, expect about 5% of 20,000 genes to be DE at fold change 1.5, and dispersion around 0.2. How many samples per group do I need?"

> "I have a 4-sample pilot. Estimate dispersions with DESeq2 and use them to size the full study."

### Single-cell

> "How should I split a fixed scRNA-seq budget between number of patients and cells per patient for differential expression between conditions?"

### Budget and Verification

> "My grant only funds 6 per group. What power and true FDR does that give for 1.5-fold changes?"

> "A reviewer asked why n=3 is not enough; give me a defensible sample-size justification."

## Related Skills

- power-analysis - The power-given-n direction and simulation-based power
- randomization-blocking - The experimental unit defines what is counted as a replicate
- batch-design - Balanced designs assume equal n per group
- differential-expression/deseq2-basics - Estimating pilot dispersions
- single-cell/preprocessing - Pseudobulk aggregation underlying scRNA-seq cohort sizing
- clinical-biostatistics/power-and-sample-size - Sample size for regulated clinical trials
