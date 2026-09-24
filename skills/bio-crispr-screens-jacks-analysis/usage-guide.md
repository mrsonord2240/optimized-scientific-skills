# JACKS Analysis - Usage Guide

## Overview

Run JACKS (Allen et al 2019 Genome Research) for joint Bayesian decomposition of CRISPR-screen log-fold-changes into per-sgRNA guide-efficacy and per-condition per-gene effect. Designed for multi-screen joint analysis where guide efficacy is shared (same library, same chemistry, cross-cell-line or cross-condition), enabling ~2.5x smaller screens (fewer replicates/guides) and ~21% lower gene-effect error vs MAGeCK on the same data. Outputs per-gene posterior effect + posterior std and per-sgRNA efficacy + std.

## Prerequisites

Install, required input files and their formats are in SKILL.md ("Version Compatibility" and "Run JACKS Joint Analysis"). Optional: a matched-library reference efficacy prior (SKILL.md, "Build Library-Wide Efficacy Prior from Reference Screens").

## Quick Start

Tell the AI agent what to do:
- "Run JACKS jointly on my 4 screens (same Brunello library, different cell lines) and output per-line gene effects + shared sgRNA efficacy"
- "Build a guide-efficacy prior from DepMap CRISPR screens and reuse it for my small custom-library screen to reduce sample-size needs"
- "Compare JACKS gene effects vs MAGeCK RRA on the same data; identify where they disagree and why"
- "Identify the low-efficacy sgRNAs (X1 <0.3) in my library for re-design in v2"
- "Decide whether JACKS or Chronos is right for my 10-cell-line cancer dependency screen"

## Example Prompts

### Multi-Screen Joint Analysis

> "I have 4 Brunello screens across HCT116, HEK293T, A375, and MCF7 cell lines, each with 3 replicate Day 0 and Day 14 samples. Run JACKS jointly with the default settings (`--apply_w_hp` off). Output per-line gene effects, the matching posterior-std file, and shared guide efficacy."

> "My screens were done across 6 weeks in two batches. Set up a per-batch JACKS run and compare to a single joint run; quantify whether batch sharing improves or degrades signal."

### Reference Efficacy Transfer

> "Extract the per-sgRNA efficacy posterior from a JACKS run on the DepMap Brunello panel (50 cell lines, ~10,000 screen days). Save as `brunello_efficacy_prior.tsv`. Then run JACKS on my single Brunello screen passing it via `--reffile`; efficacy-aware testing enables the ~2.5x smaller screens reported in Allen 2019."

### Library Calibration / Re-design

> "Output the per-sgRNA efficacy from JACKS analysis of my custom screen. Flag sgRNAs with X1 <0.3 as low-efficacy. Group by gene; flag any gene where all guides are low-efficacy for library re-design."

> "I want to design a v2 of my custom library. Use the JACKS efficacy output to drop the bottom 25% of guides and replace with new candidates from CRISPOR."

### Comparison and Diagnostics

> "Compute Spearman ρ between JACKS X1 and MAGeCK neg|lfc on the same dataset. Investigate any rank disagreement >100 positions in the top-1000 hit list."

> "Diagnose why JACKS shows median efficacy 0.18 across my whole library. Is this a chemistry mismatch (CRISPRi screen run with Cas9 defaults), an over-shrinkage from `--apply_w_hp`, or a real library quality issue?"

> "My JACKS p-values change between repeated runs even though the gene effects do not. Explain why and make them reproducible."

## What the Agent Will Do

1. Verify input file formats: sgRNA names consistent between counts, guide_map, replicate_map
2. Decide whether joint analysis is appropriate: same library, same chemistry, ≥3 screens
3. If applicable, build the reference efficacy prior from a matched public dataset
4. Run JACKS via `python run_JACKS.py` from `JACKS/jacks/`, or programmatically via `jacks.jacks_io.runJACKS`
5. Leave `--apply_w_hp` off unless deliberately using the hierarchical gene-effect prior (the tool's help advises caution)
6. Check convergence: genes ending at the 50-iteration cap in the DEBUG log get a refit with a higher `n_iter`; seed `random` and set `PYTHONHASHSEED` when p-values must be reproducible
7. Generate the gene-effect matrix plus its std file, and the sgRNA efficacy file (`sgrna`, `X1`, `X2`)
8. Call hits on effect/std (abs >2); supply `--ctrl_genes` (and `n_pseudo` > 0 in Python) if p-values are needed
9. Flag low-efficacy guides (X1 <0.3) and genes where all guides are weak (re-design candidates)
10. Cross-validate with MAGeCK / BAGEL2: identify high-confidence hits in agreement, single-tool hits flagged for orthogonal validation
11. Decide if Chronos is preferred (cancer-line multi-cell-line screens with CN bias)
12. Report the JACKS output files (`<outprefix>_gene_JACKS_results.txt`, `_gene_std_JACKS_results.txt`, `_grna_JACKS_results.txt`, `_gene_pval_JACKS_results.txt` only with `--ctrl_genes`, `_logfoldchange_means.txt`/`_logfoldchange_std.txt` unless `--reffile` was used, and the `.pickle`), plus, as agent-written summaries rather than JACKS outputs, the low-efficacy guide list for re-design and the JACKS-vs-MAGeCK comparison table

## Related Skills

- crispr-screens/mageck-analysis - MAGeCK MLE for joint multi-condition design (alternative)
- crispr-screens/bagel-essentiality - BAGEL2 for essentiality classification without efficacy
- crispr-screens/copy-number-correction - Chronos for cancer-line multi-screen + CN bias
- crispr-screens/screen-qc - Replicate Pearson + plasmid Gini gate JACKS use
- crispr-screens/library-design - Use JACKS efficacy output to refine library v2
- crispr-screens/hit-calling - Cross-method decision tree and reconciliation
