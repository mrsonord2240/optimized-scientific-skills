# Mediation Analysis - Usage Guide

## Overview

Decomposes the total effect of an exposure (genotype, treatment, environmental factor) on an outcome into direct and indirect paths through one or more mediators. Covers single-mediator observational mediation (`mediation::mediate`), 4-way decomposition with exposure-mediator interaction (CMAverse), high-dimensional EWAS / transcriptome-wide mediator screening (HIMA / HIMA2 / BAMA), MR-based mediation when sequential ignorability is implausible (two-step MR, MVMR-mediation), and doubly-robust double-ML mediation (`causalweight::medDML`). Sequential ignorability is fundamentally untestable, so every reported result is paired with a sensitivity analysis (Imai rho or mediational E-value).

## Prerequisites

Install commands and per-package compute time: SKILL.md `## Tool Install Notes`.

## Quick Start

Tell the AI agent what kind of mediation question is being asked:
- "Test whether expression of GENE_X mediates the effect of rs12345 on disease risk; include sensitivity to unmeasured confounding"
- "Run 4-way decomposition for treatment-mediator-outcome with exposure-mediator interaction"
- "Screen all 450k CpGs for mediators of smoking-lung-cancer association using HIMA2"
- "Run MR-mediation to estimate the fraction of the BMI-CHD effect that goes through LDL"
- "Use double-ML mediation to estimate ACME with random-forest nuisance models"
- "Compute the mediational E-value for an ACME of 0.04 with 95% CI [0.01, 0.07]"

## Example Prompts

### Single-Mediator Observational

> "I have individual-level data on genotype, gene expression, and binary disease status with covariates age/sex/PCs. Test whether expression mediates the genotype-disease association and report ACME with 95% CI and proportion mediated."

> "Run mediation analysis for SNP rs7412 with LDL cholesterol as the mediator and Alzheimer's disease as the outcome; include Imai sensitivity analysis."

### 4-Way Decomposition with Interaction

> "Decompose the smoking-COPD effect into CDE, PIE, INTref, INTmed using FEV1 as the mediator and CMAverse; outcome is binary with rare-disease prevalence around 5%."

> "Fit a regression-based 4-way decomposition with exposure-mediator interaction and exposure-set-to-0 reference."

### High-Dimensional EWAS / Transcriptome

> "Run HIMA2 on the Illumina EPIC methylation matrix to identify CpGs mediating the prenatal-smoking-birthweight relationship; control FDR at 0.05."

> "I have RNA-seq counts for 20k genes, n=200 subjects, and a binary case/control outcome. Use HIMA-binomial to find transcripts that mediate exposure-disease."

> "Run BAMA with Bayesian shrinkage on 5000 candidate protein mediators."

### MR-Mediation

> "Run two-step MR for BMI -> diastolic blood pressure -> coronary artery disease using independent instruments at each step; apply Steiger filter."

> "Use MVMR-mediation to estimate the direct effect of LDL on CHD adjusting for HDL; report conditional F for each exposure and the indirect path through HDL."

### Longitudinal / Time-Varying Confounding

> "Run g-formula mediation via CMAverse for treatment exposure with time-varying confounders measured at three timepoints."

### Sensitivity

> "Compute the mediational E-value for my ACME of 0.05 (95% CI 0.02-0.08) on a binary outcome."

> "Run medsens on this mediate result and tell me whether rho_crit is above 0.3."

For the analytic-regime decision procedure, bootstrap/FDR/F-statistic thresholds, and failure-mode handling the agent applies, see SKILL.md `## Decision Tree by Scenario`, `## Quantitative Thresholds`, and `## Common Errors`.

## Related Skills

- causal-genomics/mendelian-randomization - IV-based causal inference; foundation for MR-mediation
- causal-genomics/colocalization-analysis - Confirm shared causal variant before causal mediation
- causal-genomics/fine-mapping - Identify the causal variant driving the exposure
- methylation-analysis/differential-cpg-testing - Per-CpG inputs for HIMA EWAS mediation
- differential-expression/deseq2-basics - Expression inputs for eQTL mediation
- multi-omics-integration/mofa-integration - Multi-layer mediator construction
- population-genetics/association-testing - GWAS summary statistics for MR-mediation
- clinical-biostatistics/effect-measures - Risk-ratio / odds-ratio scales for binary outcomes
- machine-learning/model-validation - Cross-fitting and sample splitting for medDML
