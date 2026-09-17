# Mendelian Randomization - Usage Guide

## Overview

Mendelian randomization (MR) estimates the causal effect of an exposure on an outcome using genetic variants as instrumental variables (IVs). Because alleles are randomized at conception and fixed across the lifecourse, well-chosen IVs are independent of confounders that bias observational associations, allowing causal inference from GWAS summary statistics rather than randomized trials.

This skill covers the full sensitivity battery (IVW, MR-Egger, weighted median/mode, MR-RAPS, MR-PRESSO, CAUSE, MVMR, MR-Clust, LCV, LHC-MR), the bias structure of one-sample vs two-sample designs, drug-target / cis-MR for pQTL/eQTL exposures, and STROBE-MR-compliant reporting. The agent selects the primary method based on the experimental regime, runs a concordant sensitivity battery, and flags weak-IV / pleiotropy / sample-overlap issues automatically.

## Prerequisites

See SKILL.md's "Tool Installation Notes" for the install commands (CRAN + GitHub packages) and
"Common Errors" for the OpenGWAS JWT-token setup and local-clumping fallback.

## Quick Start

Tell the AI agent what to do:
- "Run a two-sample MR of LDL cholesterol on coronary heart disease from local GWAS files"
- "Test whether circulating IL-6R protein causally affects rheumatoid arthritis using cis-pQTL instruments"
- "Run multivariable MR of BMI on T2D adjusting for waist-hip-ratio"
- "Estimate the indirect effect of BMI on CHD mediated through LDL cholesterol"
- "Distinguish a causal BMI to CHD effect from correlated horizontal pleiotropy using CAUSE"
- "Produce a STROBE-MR-compliant sensitivity battery for my exposure-outcome pair"

## Example Prompts

### Standard Two-Sample MR

> "I have GWAS summary statistics for systolic blood pressure (exposure) and stroke (outcome) from independent cohorts. Run TwoSampleMR with IVW primary, plus MR-Egger, weighted median, weighted mode, MR-PRESSO, and Steiger directionality. Report effect sizes with 95% CIs and flag any pleiotropy concerns."

> "Test the causal effect of educational attainment on Alzheimer's disease using GWAS summary stats. Use local plink clumping and apply Steiger filtering."

### Drug-Target / cis-MR

> "Run a cis-MR analysis of PCSK9 inhibition on coronary artery disease using cis-pQTL instruments from UKB-PPP within 500 kb of the PCSK9 gene. Cross-validate with coloc PP.H4."

> "I want to predict the effect of an IL-6 receptor antagonist on rheumatoid arthritis. Use cis-pQTLs for IL6R from deCODE plus colocalization."

### MVMR and Mediation

> "Run multivariable MR of LDL, HDL, and triglycerides on CHD jointly. Report conditional F for each lipid trait. Use MVMR's `strength_mvmr()`."

> "Test whether BMI's effect on T2D is mediated through fasting glucose. Implement the two-step MR and MVMR difference methods."

### Non-Linear MR

> "Test for a J-shaped relationship between alcohol consumption and all-cause mortality using doubly-ranked stratification (Tian 2023)."

### Polygenic Exposure with CHP Concern

> "Run CAUSE on BMI -> coronary heart disease to separate a causal effect from correlated horizontal pleiotropy. Compare against IVW and MR-PRESSO."

### One-Sample / Sample-Overlap Concern

> "Both my exposure and outcome GWAS are from UK Biobank. Apply MRlap for joint correction of sample overlap, winner's curse, and weak-IV bias (Mounier 2023); fall back to MR-RAPS plus Burgess 2016 correction if LDSC h^2 < 0.05."

### Binary Outcomes

> "Run MR of LDL on T2D and explicitly handle the non-collapsibility of the logistic OR; report the per-SD genetically-predicted exposure effect on the marginal log-OR scale."

## What the Agent Will Do

Diagnoses the experimental regime (one-sample, two-sample, partial overlap, drug-target, polygenic,
non-linear) per SKILL.md's "Decision Tree by Experimental Scenario", then runs that regime's primary
method plus its sensitivity battery per "TwoSampleMR Standard Workflow" / "MR-PRESSO Outlier
Detection" / "CAUSE" / "MVMR with Conditional F" / "Bidirectional and Steiger", and reports per
"STROBE-MR Reporting".

## Related Skills

causal-genomics/colocalization-analysis - Required for cis-MR drug-target work to confirm a shared causal variant
causal-genomics/pleiotropy-detection - Deep MR-PRESSO / Egger / contamination-mixture diagnostics
causal-genomics/fine-mapping - Credible-set construction at instrument loci
causal-genomics/mediation-analysis - Two-step MR and MVMR difference method for X -> M -> Y
causal-genomics/transcriptome-wide-association - TWAS as MR-adjacent framework for gene-level inference
causal-genomics/proteome-mr-drug-target - Cis-pQTL drug-target MR workflow with UKB-PPP/deCODE/Fenland and coloc triangulation
population-genetics/association-testing - GWAS source for instrument selection
clinical-databases/clinvar-lookup - Annotate instrument SNPs for downstream interpretation
