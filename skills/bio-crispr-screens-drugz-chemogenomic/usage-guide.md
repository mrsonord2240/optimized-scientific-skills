# drugZ Chemogenomic - Usage Guide

## Overview

Analysis of CRISPR drug-modifier screens with drugZ (Colic et al. 2019 *Genome Medicine*): bidirectional Z-scores comparing drug vs vehicle (NOT Day-0) to identify sensitizing (synthetic-lethal) and resistance-conferring (suppressor) genes. Install, inputs, commands, thresholds, failure modes and the drugZ vs MAGeCK MLE choice are all in `SKILL.md`.

## Example Prompts

### Standard Drug Screen

> "Run drugZ on counts.txt with vehicle samples Veh_r1, Veh_r2, Veh_r3 (passed via -c flag) and drug samples Drug_r1, Drug_r2, Drug_r3 (via -x flag). Pseudocount 5 (-p 5). Output ranked by fdr_synth (sensitizers) and fdr_supp (suppressors)."

> "Run drugZ on a PARPi (olaparib) chemogenomic screen: vehicle vs olaparib at 14 days. Expect sensitizers in DDR pathway (BRCA1, BRCA2, RAD51, FANCD2), PARP1 as the expected suppressor. Identify novel resistance genes."

### Reference Choice

> "Diagnose why my drugZ output has no sensitizers. Confirm: am I comparing drug vs vehicle, or drug vs Day 0? If Day-0, re-run with vehicle samples."

### Dose Response

> "For my 3-dose PARPi screen, run drugZ at low, mid, and high dose vs vehicle. Identify dose-consistent sensitizers (same gene appearing at all 3 doses with same direction)."

### Removing Essentials

> "Run drugZ with `-r` set to a comma-delimited list of CEGv2 essential gene names. Compare to default output; the difference is the essential genes removed from the analysis."

### Method Comparison

> "Compare drugZ vs MAGeCK MLE on the same drug screen. Compute Jaccard similarity at FDR 0.05. For drugZ-only hits, investigate per-sgRNA evidence."

> "MAGeCK and drugZ disagree on PARPi resistance hits. drugZ shows 12 hits at fdr_supp <0.05; MAGeCK shows 4 at pos|fdr <0.05. Reconcile."

### Diagnostics

> "drugZ hits are dominated by RPS / RPL / EIF essentials. Diagnose: am I running with Day-0 baseline by mistake, or do I need to exclude essentials?"

> "Replicate-to-replicate drugZ runs give different top hits. Diagnose: insufficient sgRNAs per gene, or biological variation between replicates?"

## Related Skills

- crispr-screens/mageck-analysis - MAGeCK MLE alternative
- crispr-screens/bagel-essentiality - Tumor-suppressor-sensitive alternative
- crispr-screens/hit-calling - Cross-method decision tree
- crispr-screens/screen-qc - Pre-drugZ replicate concordance
- crispr-screens/library-design - 6+ sgRNAs/gene library
- crispr-screens/copy-number-correction - Pre-correction for cancer-line drug screens
- crispr-screens/base-editing-analysis - Variant-function drug screens
- pathway-analysis/go-enrichment - Functional analysis of drug-modifier hits
