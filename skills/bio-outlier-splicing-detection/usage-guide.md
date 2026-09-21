# Outlier Splicing Detection - Usage Guide

## Overview
Detect aberrant splicing in single rare-disease patients vs a control panel (research use, not a clinical report). The question is fundamentally different from differential splicing between groups. Tools: FRASER 2 (splicing outliers), OUTRIDER (expression outliers), LeafcutterMD (annotation-free intron outliers) and DROP (Snakemake pipeline integrating FRASER2 + OUTRIDER + monoallelic expression). Commands, cohort-size numbers, thresholds and licences are in `SKILL.md`.

## Prerequisites
Install notes are in `SKILL.md` (Install).

## Example Prompts

### Single-Patient Workflow
> "I have RNA-seq from a rare-disease patient and 50 control samples; run FRASER 2 with Intron Jaccard Index and report aberrant junctions with padj<0.05 and |delta|>=0.1."

### Pipeline Setup
> "Configure DROP for our research cohort with the patient + 80 in-house controls."

### Variant Follow-up
> "I have a SpliceAI hit at chr5:1234567C>T in patient X; check whether FRASER2 detects an aberrant junction within 1kb in the same sample."

### Disease-Specific
> "For ALS post-mortem brain RNA-seq, use FRASER2 to detect TDP-43-loss cryptic exons (UNC13A, STMN2, ATG4B)."

### Hyperparameter Tuning
> "Run estimateBestQ on my cohort to determine the latent dimension q for FRASER2."

### Expression and Annotation-Free Outliers
> "Use OUTRIDER to find genes with aberrant expression in the patient" / "Run LeafcutterMD for annotation-free outlier intron usage."

## Related Skills

- splice-variant-prediction - SpliceAI / Pangolin for in-silico prediction first
- differential-splicing - When testing multiple patients vs controls
- splicing-qc - Library / depth / tissue prerequisites
- variant-calling/clinical-interpretation - ACMG/AMP framework (clinical-laboratory step, outside this Skill)
- workflows/clinical-trial-pipeline - Trial-grade RNA-seq pipelines
