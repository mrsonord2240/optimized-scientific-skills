# Base Editing Analysis - Usage Guide

## Overview

Decision-grade analysis of base-editor variant-function screens. Covers base-editor library design, the Hanna 2021 BRCA1/2 SNV-scanning methodology (*Cell* 184:1064), Cuella-Martin 2021 complementary DDR-gene saturation screen (*Cell* 184:1081-1097), CBE vs ABE chemistry selection (BE3/BE4 vs ABE7.10/ABE8.20/ABE8e), editing window math (positions 4-8 from PAM-distal end), bystander attribution strategies, editing-efficiency filtering before hit calling, indel-byproduct interpretation, and the Broad be-validation-pipeline notebooks for CRISPResso2 post-processing.

## Quick Start

Tell the AI agent what to do:
- "Design a CBE library tiling BRCA1 exons 1-23 with 10-15 sgRNAs per amino acid; flag intended target Cs and bystander Cs in the editing window"
- "Run editing-efficiency filtering on pilot library counts; drop sgRNAs <30% target editing"
- "Score variant function per amino acid; deconvolute bystander confounding via allele-frequency tables"
- "Reconcile a screen hit between Hanna's CBE methodology and orthogonal prime-editor follow-up"
- "Choose BE vs PE for installing a C->T variant: BE3 if target at pos 5 with no bystanders; PE2 if zero-bystander requirement"
- "Diagnose: my BE sample has 60% editing but 35% are indels -- is this Cas9 contamination?"

## Example Prompts

### Library Design

> "Design a CBE saturation library tiling BRCA1 RING domain (amino acids 1-100). 10-15 sgRNAs per amino acid where at least one C in the editing window (positions 4-8) hits the target codon. Annotate each sgRNA with target + bystander editable-base positions. Output library.tsv with sgRNA, target_aa, target_positions, bystander_positions columns."

> "I need to install MLH1 c.677A>G as a single intended variant. CBE won't work (need ABE). Find ABE7.10 or ABE8e sgRNAs that place A at position 5 with no bystanders. If no zero-bystander spacer exists, list candidates sorted by bystander_count and recommend prime editor as alternative."

### Editing Efficiency Filtering

> "Run CRISPResso2 on my pilot timepoint samples. Compute target editing % per sgRNA. Keep sgRNAs >30% target editing for the primary screen; flag >50% as validation-grade. Output filtered library and editing efficiency report."

> "My library shows median 22% editing across guides. Diagnose: cell-line BE activity issue, vector mismatch, or BE chemistry choice (CBE vs ABE)?"

### Hit Calling for Variant Function

> "Apply MAGeCK MLE to the BE screen counts (vehicle vs PARPi drug). Use only editing-efficient sgRNAs. Aggregate per-sgRNA LFC to per-variant scores; deconvolute bystanders via allele tables; output per-variant fitness with confidence based on the number of sgRNAs hitting each variant."

> "Compare drugZ vs MAGeCK MLE on the same BE drug-modifier screen; identify resistance and sensitivity variants using the Hanna 2021 CBE variant-scanning approach."

### Bystander Deconvolution

> "From CRISPResso2 allele tables, separate reads by edit pattern: target only, target+bystander_1, target+bystander_2. Compute per-pattern fitness contribution; flag variants where bystander dominates."

> "For BRCA1 R71 -> R71X variant, I have 5 sgRNAs with different bystander patterns. Identify variants attributable to R71X alone by finding consistent signal across diverse bystander backgrounds."

### Method Comparison

> "Compare CBE BE3 vs BE4max vs eA3A-BE3 vs evoCDA-BE for tiling BRCA1. Score: editing window width, bystander rate, indel byproduct, target-base preference (TC vs other)."

> "Cross-validate Hanna 2021 BRCA1 BE screen hits with parallel PE screen at the same variants. Report concordant variants and BE-only (likely bystander-confounded) hits."

### Diagnostics

> "My BE sample has 70% editing but 25% indels. Compute substitution-vs-indel ratio; if <3, diagnose as Cas9 contamination."

> "Allele table shows 35% target+bystander double-edit. Deconvolute: how much is target-attributable vs bystander-driven?"

## Related Skills

- crispr-screens/crispresso-editing - CRISPResso2 modes and allele tables
- crispr-screens/library-design - base-editor library design
- crispr-screens/prime-editing-screens - Orthogonal PE for variant attribution
- crispr-screens/hit-calling - Variant-level hit aggregation
- crispr-screens/screen-qc - Editing-efficiency QC
- crispr-screens/drugz-chemogenomic - drugZ for BE drug-modifier screens
- crispr-screens/mageck-analysis - MAGeCK MLE for BE screen sgRNA-level analysis
- clinical-databases/clinvar-lookup - Pathogenicity annotation
- variant-calling/variant-annotation - VEP for predicted protein effects
