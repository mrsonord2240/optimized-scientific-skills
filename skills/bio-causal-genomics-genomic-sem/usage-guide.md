# Genomic SEM Usage Guide

## Overview

Genomic SEM fits structural equation models to GWAS summary statistics, treating each GWAS as an indicator of one or more latent genetic factors. The framework (Grotzinger 2019 Nat Hum Behav 3:513) uses an LDSC-derived genetic covariance matrix (S) and its sampling covariance (V) as inputs to a `lavaan`-based SEM, supporting:

- Common-factor confirmatory models across correlated traits
- User-specified factor structures with arbitrary lavaan syntax
- Exploratory structural equation modeling (ESEM)
- Multivariate GWAS in which a SNP is regressed on a latent factor (`commonfactorGWAS`) with the Q_SNP heterogeneity test
- Stratified GenomicSEM partitioning factor heritability across functional annotations
- Cross-checks against MTAG multi-trait analysis

The skill emphasizes when GenomicSEM is the right multivariate framework, when MTAG is a better fit, and how to read fit indices and Q_SNP to avoid over-claiming "common-factor SNPs" that are in fact trait-specific.

## Prerequisites

Package versions (GenomicSEM must be paired with lavaan 0.6.19, not a newer lavaan -- this
breaks 3 of 4 core functions), install commands, and reference-data requirements: SKILL.md
`## Version Compatibility` and `## Tool Installation`.

## Quick Start

Tell an AI agent what to model:

- "Fit a common-factor model to my 5 psychiatric GWAS using GenomicSEM and report CFI and RMSEA"
- "Run a common-factor GWAS across MDD, anxiety, and PTSD; flag Q_SNP-significant SNPs"
- "Compare GenomicSEM common-factor results against MTAG on the same input traits"
- "Fit a two-factor user model: F1 loads on lipid traits, F2 loads on glycemic traits, with correlated factors"
- "Diagnose this Heywood case (negative residual variance) in my GenomicSEM output"
- "Partition the heritability of the common factor across baseline-LD annotations"
- "Check whether the sampling covariance V is positive definite before fitting"

## Example Prompts

### Common-Factor Model

> "Take these three educational-attainment GWAS sumstats, run `GenomicSEM::ldsc()` against the EUR reference, fit a common-factor CFA, and report CFI, RMSEA, SRMR, and standardized loadings. Diagnose any Heywood case if present."

### Common-Factor GWAS with Q_SNP

> "Run a common-factor GWAS across MDD, BIP, and SCZ using `commonfactorGWAS` with `DWLS` estimation. Report SNPs with factor p < 5e-8 AND Q_SNP p > Bonferroni-corrected threshold. Exclude Q_SNP-significant SNPs from the common-factor SNP list."

### Confirmatory Two-Factor Model

> "Fit a two-factor user model where F1 = LDL + HDL + triglycerides and F2 = fasting glucose + HbA1c + 2hr glucose, with F1 ~~ F2 free to estimate factor correlation. Use `usermodel` with DWLS estimation. Report fit indices and factor correlation."

### MTAG Comparison

> "Run MTAG via CLI on the same three lipid GWAS, then compare per-trait top hits with GenomicSEM common-factor SNPs. Check that MTAG `maxFDR` is < 5% for each trait. Where they disagree, classify the SNP as factor-mediated vs trait-specific."

### Stratified GenomicSEM

> "Use `s_ldsc()` to partition the heritability of the latent externalizing factor across baseline-LD annotations + ENCODE/Roadmap cell-type marks. Report per-annotation factor enrichment and the cell type with the largest factor-tau coefficient."

### Sample Overlap Diagnosis

> "These four GWAS are from UK Biobank. Verify the bivariate LDSC intercepts and confirm the V matrix off-diagonals reflect overlap. Fit the common-factor model and check that SEs differ from a naive analysis that ignores V."

### ESEM Exploratory

> "Run an ESEM with two free factors on these 8 cognition + personality GWAS to discover whether a single g factor captures the shared variance or whether g + p_factor emerges. Report rotation: geomin oblique."

For the standard end-to-end procedure (munge -> ldsc -> fit -> inspect -> Q_SNP -> MTAG ->
reconcile), fit-index thresholds, Heywood/overlap/Q_SNP/MaxFDR failure-mode handling, and
identification rules the agent applies, see SKILL.md `## Standard Workflow`, `## Model Fit
Diagnostics`, `## Quantitative Thresholds`, `## Per-Method Failure Modes`, and `## Common
Errors`.

## Related Skills

causal-genomics/mendelian-randomization - Run MR with the common-factor sumstats as exposure
causal-genomics/heritability-partitioning - sLDSC foundations + LDAK comparison; required upstream for stratified GenomicSEM
causal-genomics/colocalization-analysis - Resolve overlap of factor-significant SNPs with eQTL/pQTL signals
causal-genomics/pleiotropy-detection - Q_SNP is the per-SNP pleiotropy diagnostic in GenomicSEM
causal-genomics/fine-mapping - Construct credible sets at factor-significant loci
causal-genomics/mediation-analysis - SEM mediation overlaps with `usermodel` indirect path coefficients
causal-genomics/transcriptome-wide-association - TWAS on factor GWAS output
population-genetics/association-testing - GWAS workflow upstream of GenomicSEM
