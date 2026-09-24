# Batch Correction - Usage Guide

## Overview

Decision-grade batch correction for pooled CRISPR screens. Covers diagnosis (PCA + variance decomposition to decide whether batch dominates), the four primary correction strategies (ComBat empirical-Bayes, RUV, SVA, NTC-anchored normalization), the model-based alternative of including batch as a MAGeCK MLE / Chronos covariate (preferred for hit calling), screen-specific batch sources (cell-passage cohort, library lot, infection day, sequencing run, FBS lot, Cas9 lot), and the failure modes when correction destroys biology.

## Prerequisites

Install commands, required inputs and the versions the code was checked on are in `SKILL.md` ("Version Compatibility").

## Quick Start

Tell the AI agent what to do:
- "Diagnose batch effects in my screen: PCA, variance decomposition, is correction needed?"
- "Apply ComBat correction with condition as the `mod` covariate; verify PR-AUC is preserved post-correction"
- "Add batch as a covariate to MAGeCK MLE design matrix instead of pre-correcting"
- "Compare ComBat vs NTC-anchored normalization on my multi-batch screen"
- "Decide if I should correct or redesign because batch is fully confounded with condition"

## Example Prompts

### Diagnostics

> "Run PCA on my multi-batch screen colored by batch and condition. Compute F-statistics for batch and condition in each PC. Decide whether batch dominates (correction needed) or condition dominates (no correction needed)."

> "Replicate Pearson within batches is 0.95+ but across-batch is 0.78. Confirm batch effects are present; recommend ComBat vs NTC-anchored vs covariate modeling."

> "Check whether batch is confounded with condition. If yes (e.g., all drug arm in batch 2), correction will destroy biology -- explicit covariate modeling required instead."

### ComBat Application

> "Apply ComBat to log2(counts+1) with condition as the `mod` covariate. Verify post-correction PCA shows batches overlap and CEGv2 PR-AUC is preserved."

> "Compare ComBat with mod (biological covariate) vs without. Quantify how much biological signal is preserved in each case."

### Batch-Aware MLE

> "Build a MAGeCK MLE design matrix with explicit batch indicator columns. Run mageck mle on the multi-batch screen and report per-condition beta scores after batch adjustment. Compare to ComBat-then-MAGeCK approach."

> "Add SVA-discovered latent factors as covariates in MAGeCK MLE design matrix to capture unknown batch confounders."

### NTC-Anchored Normalization

> "Use the 800 non-targeting controls in my library as normalization anchor. Scale each sample so NTC median = 1000. Then run hit calling on NTC-anchored counts."

### Multi-Cell-Line Cancer Panels

> "I have 12 cancer cell lines screened across 4 batches. Switch to Chronos (built-in batch + CN modeling) instead of ComBat + MAGeCK; explain why this is preferred."

### Failure Diagnostics

> "After ComBat, my essentialome PR-AUC dropped from 0.78 to 0.41. Diagnose: did ComBat eliminate biological signal? Re-run with mod covariate or revert and use MLE-with-batch-covariate."

## What the Agent Will Do

Diagnose (PCA + variance decomposition), check batch-versus-condition confounding, choose the correction, re-check PCA and CEGv2 PR-AUC afterwards. The decision tree, thresholds, failure modes and the validation checklist are all in `SKILL.md`.

## Related Skills

- crispr-screens/mageck-analysis - Batch-aware MLE design matrix (preferred)
- crispr-screens/screen-qc - Pre-correction PCA + variance decomposition
- crispr-screens/copy-number-correction - Chronos handles batch + CN jointly
- crispr-screens/library-design - NTC composition for NTC-anchored norm
- crispr-screens/hit-calling - Post-correction hit calling
- crispr-screens/jacks-analysis - Joint analysis across batches
- crispr-screens/in-vivo-screens - In-vivo-specific batch sources
