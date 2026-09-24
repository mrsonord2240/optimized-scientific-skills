# BAGEL2 Essentiality - Usage Guide

## Overview

Decision-grade essentiality calling for CRISPR-Cas9 fitness screens using BAGEL2 (Kim & Hart 2021 *Genome Medicine*). Computes Bayes Factors per gene by log-likelihood ratios over per-sgRNA fold changes, anchored against CEGv2 reference core-essentials (684 genes) and NEGv1 reference non-essentials (927 genes). Improvements over BAGEL1: linear extrapolation enables tumor-suppressor detection; multi-target off-target correction reduces false positives. Thresholds are in SKILL.md.

## Prerequisites

BAGEL2 from `git clone https://github.com/hart-lab/bagel` (no PyPI release), the CEGv2/NEGv1 reference files from that repo, and a count matrix (sgRNA, GENE, sample columns) with a Day 0 or plasmid control column.

## Quick Start

Tell the AI agent what to do:
- "Run BAGEL2 on my Brunello screen: compute fold changes then Bayes factors using CEGv2/NEGv1; pick BF threshold from precision-recall analysis"
- "Identify tumor-suppressor candidates from BAGEL2 negative BF values; cross-check against COSMIC tumor suppressor list"
- "Compare BAGEL2 BF vs MAGeCK FDR on the same dataset; reconcile disagreements at BF 5-7 boundary"
- "Diagnose BAGEL2 calling no hits despite known essentials -- is it a reference-set issue?"
- "Calibrate a high-stringency BF threshold: BF >12 vs BF >30, using BAGEL.py pr on my own screen"

## Example Prompts

### Standard Workflow

> "Run BAGEL.py fc then bf on counts.txt with controls Plasmid and treatment samples Drug_r1,Drug_r2,Drug_r3. Output bayes_factor.txt sorted by BF descending."

> "Run BAGEL.py pr after to generate precision-recall curve against CEGv2; pick the BF threshold that reaches 95% precision on my screen's own PR curve."

### Tumor Suppressor Detection

> "From BAGEL2 output, identify genes with BF <-6 (significantly negative). Cross-reference against COSMIC tumor suppressor gene list. Output candidate tumor suppressors with their BF, fold change, and sgRNAs."

> "BAGEL2 calls 80 tumor-suppressor candidates in my dropout screen. Many don't replicate in literature -- investigate whether my screen design supports tumor-suppressor calling (drug-modifier vs simple dropout)."

### Calibration and Comparison

> "Compare BAGEL2 BF >6 hits vs MAGeCK RRA neg|fdr <0.05 hits on the same screen. Compute Jaccard similarity. Investigate top 10 disagreements."

> "Adjust BF threshold based on screen quality. High-quality screen (PR-AUC against CEGv2 >0.85) supports BF >6. Lower-quality may need BF >12 for same precision."

### Diagnostics

> "My BAGEL2 returns no genes with BF >0. Diagnose: wrong reference set file, library coverage too low, or genuinely flat screen?"

> "Per-sgRNA LLR contributions show one guide dominating BF for several hits. Apply second-best-sgRNA rule from [[hit-calling]] to filter these guide-of-one hits."

## More

Commands, thresholds, seed handling, failure modes and interpretation rules live in SKILL.md (the pre-flight and post-run checks are `examples/check_bagel_inputs.py`; the end-to-end script is `examples/run_bagel2.sh`).

## Related Skills

- crispr-screens/mageck-analysis - MAGeCK as alternative
- crispr-screens/drugz-chemogenomic - drugZ for drug screens (also tumor-suppressor sensitive)
- crispr-screens/jacks-analysis - JACKS efficacy diagnostics for guide-of-one
- crispr-screens/hit-calling - Cross-method consensus
- crispr-screens/screen-qc - Pre-BAGEL CEGv2 PR-AUC
- crispr-screens/library-design - 4-6 sgRNAs/gene library standard
- crispr-screens/copy-number-correction - Pre-correction for cancer-line screens
- pathway-analysis/go-enrichment - Downstream
