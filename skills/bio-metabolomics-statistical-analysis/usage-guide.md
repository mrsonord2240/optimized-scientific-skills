# Statistical Analysis Usage Guide

## Overview

Statistical analysis decides which metabolites differ between conditions and whether a discriminant model is real. The two endemic failure modes it guards against are overfitting (in p >> n, supervised models separate pure noise -- a score plot proves nothing without permutation-validated Q2) and false-discovery inflation (thousands of correlated tests on small n yield long "significant" lists that do not replicate). It covers transformation/scaling as an explicit modeling choice, unsupervised structure (PCA/HCA for QC), validated PLS-DA/OPLS-DA, univariate testing with covariate adjustment, and dependence-aware multiple testing.

## Prerequisites

Install commands and the conceptual checks to settle first are in SKILL.md (top of the file).

## Quick Start

Tell your AI agent what you want to do:
- "Run PCA with Pareto scaling and check whether my pooled QC samples cluster tightly"
- "Fit an OPLS-DA model and permutation-test it -- I only trust it if pQ2 is small"
- "Run Welch t-tests between case and control with explicit BH FDR and fold changes"
- "Fit a linear mixed model per metabolite for my longitudinal design"
- "Show me whether my top hits survive switching from Pareto to unit-variance scaling"

## Example Prompts

### Scaling and Unsupervised Structure
> "Transform and Pareto-scale my feature table, run PCA, and color the scores by batch and injection order to check for drift."
> "Re-run the PCA with unit-variance scaling and tell me whether the top loadings change."

### Validated Multivariate Models
> "Build an OPLS-DA model between disease and control, run 1000 permutations, and report R2X, R2Y, Q2, and pQ2."
> "Put a PCA score plot next to my PLS-DA plot so I can see whether the separation is real or supervised-only."
> "Extract VIP scores but only after the model passes the permutation test, and cross-check them against my univariate hits."

### Univariate Testing and FDR
> "Run Welch t-tests in Python with Benjamini-Hochberg correction and report log2 fold changes."
> "Fit a per-metabolite linear model adjusting for age, sex, and BMI."
> "My features are heavily correlated -- use an effective-number-of-tests correction instead of Bonferroni."

### Reproducibility
> "Check whether my discriminant result survives changing the scaling and the imputation method."
> "I have a validation cohort -- estimate external performance, not just cross-validated discovery performance."

## Related Skills

- metabolomics/normalization-qc - Sample-wise normalization, drift/batch correction, missing-value imputation upstream of testing
- metabolomics/pathway-mapping - Functional interpretation of differential metabolites
- machine-learning/biomarker-discovery - Feature selection inside CV, stability, minimal-optimal vs all-relevant
- machine-learning/model-validation - Leakage taxonomy, nested CV, calibration vs discrimination
- experimental-design/multiple-testing - FDR vs FWER regime, discovery vs confirmatory
- data-visualization/volcano-and-ma-plots - Volcano plot recipes

References are listed in SKILL.md.
