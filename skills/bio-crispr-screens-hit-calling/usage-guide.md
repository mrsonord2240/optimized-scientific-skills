# Hit Calling - Usage Guide

## Overview

Decision-grade cross-method orchestration for calling significant hits in pooled CRISPR screens. Catalogs the seven primary methods (MAGeCK RRA, MAGeCK MLE, BAGEL2, drugZ, JACKS, Chronos, CERES), the experimental design each is built for, the failure modes outside their domain, and the reconciliation logic when methods disagree. Defines tier-based confidence stratification (Tier 1 = 3-method consensus; Tier 2 = 2-method; Tier 3 = single-method exploratory).

## Prerequisites

```bash
# All five primary tools
conda install -c bioconda mageck                     # RRA + MLE (not on PyPI)
git clone https://github.com/hart-lab/bagel          # BAGEL2 (PyPI 'bagel' is an unrelated package)
git clone https://github.com/felicityallen/JACKS     # JACKS (PyPI 'jacks' is an unrelated package)
git clone https://github.com/hart-lab/drugz          # drugZ (not on PyPI)
# Chronos (DepMap)
pip install crispr_chronos
# Custom analyses
pip install pandas numpy scipy statsmodels matplotlib seaborn
```

## Quick Start

Tell the AI agent what to call:
- "Decide which hit-calling method to run on my screen design (2-condition vs time-course vs cancer-line vs drug screen)"
- "Run MAGeCK + BAGEL2 + drugZ on the same data and build a tier-1/tier-2/tier-3 consensus hit list"
- "Reconcile MAGeCK and BAGEL2 disagreement on my essentiality screen"
- "Apply the second-best-sgRNA rule to my hit list to flag single-guide-driven false positives"
- "Pick FDR threshold by screen quality: high-quality (PR-AUC >0.85) -> FDR 0.05; low-quality (PR-AUC 0.5-0.7) -> FDR 0.01"

## Example Prompts

### Method Selection by Design

> "My screen is 5 cancer cell lines vs Day 0 controls across 14 days. Pick Chronos vs MAGeCK MLE vs JACKS and explain why."

> "I have one drug screen: vehicle vs drug at MOI 0.3, 14 days, 3 replicates each. Pick drugZ vs MAGeCK MLE."

> "I have a time-course screen with Day 0 / 7 / 14 / 21. Should I use MAGeCK MLE or run RRA pairwise?"

> "My screen is a combinatorial paired-guide library testing GIs. Decide on MAGeCK MLE with GI scoring vs custom analysis."

### Multi-Method Consensus

> "Run MAGeCK RRA and BAGEL2 on my essentiality screen. Build a 2-method consensus hit list. Use BF >6 as BAGEL2 threshold (≈ 90% posterior; ~5% FDR by convention). Output Tier 1 (both methods) and Tier 2 (one method) lists separately."

> "Run MAGeCK + BAGEL2 + drugZ on my drug screen. Output the tier-1 consensus (3-method agreement) at FDR <0.05 / BF >6 across all. These hits go to arrayed validation."

> "Compute Spearman ρ between MAGeCK neg|score and BAGEL2 BF, sign-correcting first (`-neg|score` vs `BF` -- see SKILL.md's Correlating MAGeCK and BAGEL2 Scores section). If the sign-corrected ρ <0.6, audit why the two methods disagree." (Naive, uncorrected ρ on real data is -0.81 and would falsely read as disagreement; sign-corrected ρ is +0.81 -- strong agreement.)

### Reconciliation

> "MAGeCK says ERBB2 is essential in HER2+ SK-BR-3 cells. Chronos says it isn't. Reconcile -- is this real essentiality or copy-number artifact?"

> "JACKS calls 200 hits that MAGeCK does not. Sample the top 20 disagreements and explain which to trust based on per-sgRNA efficacy."

> "BAGEL2 returns BF >12 for a gene that MAGeCK gives neg|fdr >0.5. Investigate the per-sgRNA fold changes; identify which method is correct."

### Threshold Selection

> "Suggest the right hit threshold for my screen quality grade. Pass-quality (CEGv2 PR-AUC >0.85) -> standard FDR; CN-confounded -> apply Chronos; low-quality -> use BAGEL2 with calibrated BF instead of FDR."

> "My screen has PR-AUC 0.72 (passing but not excellent). Should I tighten FDR from 0.05 to 0.01? Make the tradeoff explicit between sensitivity and specificity."

### Custom Hit Calling

> "Run a z-score hit calling using only non-targeting controls as the null. Output gene-level z, p-value, BH-adjusted FDR."

> "Apply the second-best-sgRNA rule: a gene is a hit only if the 2nd-most-extreme sgRNA also passes the threshold. Filter MAGeCK output by this rule."

## Where the method lives

The decision tree, thresholds, reconciliation rules, confidence tiers, failure modes and the order of operations are all in SKILL.md; this guide does not repeat them.

## Related Skills

- crispr-screens/mageck-analysis - MAGeCK RRA + MLE details
- crispr-screens/bagel-essentiality - BAGEL2 details
- crispr-screens/drugz-chemogenomic - drugZ for drug screens
- crispr-screens/jacks-analysis - JACKS for multi-screen
- crispr-screens/copy-number-correction - Chronos / CRISPRcleanR for cancer-line CN bias
- crispr-screens/screen-qc - Quality gate before hit calling
- crispr-screens/combinatorial-screens - GI scoring
- crispr-screens/perturb-seq-analysis - SCEPTRE for single-cell
- pathway-analysis/go-enrichment - Functional analysis of hit lists
- pathway-analysis/gsea - GSEA on ranked gene lists
