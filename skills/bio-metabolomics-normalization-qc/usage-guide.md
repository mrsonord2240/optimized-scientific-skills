# Normalization and QC Usage Guide

## Overview

Untargeted-metabolomics preprocessing is a chain of assumption-laden modeling steps, not neutral cleanup: drift correction, batch alignment, sample normalization, and imputation each impose a hypothesis about where the unwanted variance lives, and a wrong hypothesis produces confidently wrong results rather than merely noisier ones. This skill guides QC design, feature filtering by D-ratio and RSD, within-batch drift correction, QC-anchored batch correction, sample normalization (PQN/MSTUS), and mechanism-aware missing-value imputation. It guards against the classic traps: over-correction that games the QC-clustering metric, TIC closure that smears one feature's change across all others, ComBat fabricating signal under imbalance, and MAR imputation erasing left-censored on/off biology. Install notes, required inputs, pipeline order, thresholds and guardrails are in SKILL.md.

## Quick Start

Tell your AI agent what you want to do:
- "Filter features by D-ratio and QC RSD, then report how many each filter removed"
- "Correct injection-order drift with QC-RSC and check it on held-out QCs"
- "PQN-normalize my urine samples and check the dilution factor isn't correlated with group"
- "Impute missing values by mechanism -- QRILC for left-censored, RF for sporadic"
- "Tell me whether ComBat is safe for my batch design or if I should use QC-anchored alignment"

## Example Prompts

### QC Design and Feature Filtering
> "Exclude the conditioning injections, then filter features by blank ratio and detection rate."
> "Compute the robust D-ratio for every feature and keep those below 0.5, reporting the count removed."
> "Lead my filtering with D-ratio instead of a blind 30% CV cutoff and explain why."

### Drift and Batch Correction
> "Apply QC-RSC drift correction per feature against injection order using my pooled QCs."
> "My cohort is over 800 samples with complex batch structure -- should I use SERRF instead of LOESS?"
> "Check whether ComBat is safe given my batch design, or whether biology is confounded with batch."
> "Validate the drift correction on held-out QCs and dilution-series linearity, not QC clustering."

### Sample Normalization
> "PQN-normalize these urine samples and check the dilution factor doesn't track the phenotype."
> "Compare TIC vs PQN and tell me whether closure is inflating a coordinated change."

### Missing Values
> "Diagnose whether missingness is left-censored or random, then impute accordingly."
> "Filter features missing in more than half the samples before imputing the residual holes."
> "Re-run my key result under QRILC and missForest to check the hits aren't an imputation artifact."

## Related Skills

- metabolomics/xcms-preprocessing - Generates the feature table this skill consumes
- metabolomics/msdial-preprocessing - Alternative feature-table source
- metabolomics/statistical-analysis - Transformation/scaling and downstream multivariate stats
- experimental-design/batch-design - Randomization and design that make correction valid
- differential-expression/batch-correction - ComBat/SVA mechanics shared with transcriptomics
