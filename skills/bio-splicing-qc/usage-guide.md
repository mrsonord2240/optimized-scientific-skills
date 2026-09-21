# Splicing QC - Usage Guide

## Overview
Checks whether RNA-seq data can support alternative splicing analysis: design audit, STAR 2-pass alignment, junction saturation, known-vs-novel junctions, junction support, splice-site strength, strandedness, 3' bias and rRNA contamination. Failures in these layers silently bias PSI estimates and inflate novel-junction false positives. Commands, thresholds and failure modes are in `SKILL.md`; `examples/splicing_qc.py` is the runnable helper and `examples/test_splicing_qc.py` its self-test.

Research QC only: scores describe a dataset, not an individual patient's variant.

## Example Prompts

### Pre-Sequencing Design Review
> "Will PE 100nt poly(A)-selected libraries at 30M reads/sample support splicing analysis? If not, what should change?"

### Post-Alignment QC
> "Run junction saturation and junction annotation across all my BAMs and report samples below acceptable thresholds."

### STAR 2-Pass
> "Configure cohort-style STAR 2-pass alignment for 12 samples - collect SJ.out.tab from pass 1 and re-align everything in pass 2 with the augmented junction set."

### Splice-Site Scoring
> "For my candidate cryptic splice sites, compute MaxEntScan donor and acceptor scores plus SpliceAI in-vivo usage probability."

### Diagnostic
> "I'm seeing a high novel-junction rate (>40%); diagnose whether it's annotation gaps, mapping artifacts, contamination, or biology (TDP-43, SF3B1)."

### Library Checks
> "Is my library stranded, rRNA-contaminated or 3' biased, and what should I pass to rMATS for strand?"

## Related Skills

- splicing-quantification - PSI estimation after QC passes
- read-alignment/star-alignment - STAR 2-pass detail
- read-qc/quality-reports - General sequencing QC
- read-qc/contamination-screening - rRNA / contamination
- splice-variant-prediction - SpliceAI / Pangolin for variant impact
- long-read-splicing - When short-read QC is fundamentally limiting
