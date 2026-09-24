# Fine-Mapping - Usage Guide

## Overview

Resolve GWAS lead SNPs to credible sets of likely causal variants by fitting sparse Bayesian regressions that propagate linkage disequilibrium (LD) into posterior inclusion probabilities (PIPs). Modern fine-mapping is dominated by SuSiE (Wang 2020) and its variants: susie_rss for summary statistics, SuSiE-inf for non-sparse loci, SuSiEx for cross-ancestry, and coloc.susie for downstream colocalization. FINEMAP, CAVIAR, DAP-G, PAINTOR, PolyFun, MultiSuSiE, and FOCUS cover specialized scenarios. The hardest practical problem is LD reference mismatch; this skill bakes in the `estimate_s_rss` and `kriging_rss` diagnostics as a mandatory step.

Install commands for every tool are in `SKILL.md` (Tool Install Notes).

## Quick Start

Tell the AI agent what to do in natural language:
- "Fine-map this GWAS locus to a 95 percent credible set using susie_rss"
- "Run the LD diagnostic estimate_s_rss before reporting credible sets"
- "Compare SuSiE and FINEMAP at this locus and reconcile disagreements"
- "Apply PolyFun functional priors to sharpen PIPs"
- "Cross-ancestry fine-map using SuSiEx with EUR, EAS, and AFR summary statistics"
- "Fine-map the HLA region with L=30 and explain why credible sets stay wide"
- "Feed susie_rss credible sets into coloc.susie for two-trait colocalization"

## Example Prompts

### Single-Locus EUR GWAS
> "Fine-map a 1 Mb window around rs12345 on chr6 using susie_rss with 1000G EUR LD. Run estimate_s_rss and report lambda. Extract 95 percent credible sets, purity, and top PIP variants."

> "Same locus, but try L=10 and L=20 and report whether the larger L changes the credible-set count."

### Polygenic / Non-Sparse Locus
> "This UK Biobank locus has 200 SNPs with -log10(p) > 4. Standard SuSiE produces 8 small credible sets that do not replicate. Refit with SuSiE-inf and compare."

### Cross-Ancestry
> "Run SuSiEx on EUR (N=500k), EAS (N=200k), and AFR (N=80k) summary statistics at chr1:50-51Mb. Compare credible set size to single-ancestry susie_rss in EUR."

### Functional Priors
> "Compute PolyFun per-SNP priors genome-wide using the baseline-LF reference, then fine-map this locus with susie_rss prior_weights set from PolyFun output. Compare to uniform-prior PIPs."

### TWAS Fine-Mapping
> "I have FUSION TWAS Z-scores for 20 co-regulated genes at chr19:45Mb. Run FOCUS to identify the likely causal gene."

### LD Diagnostic
> "Run estimate_s_rss and kriging_rss on this locus; flag any SNPs with |z_obs - z_exp| > 3 and explain whether the LD reference is suitable."

### HLA Caveat
> "I tried to fine-map a chr6:30-33 Mb autoimmune locus with susie_rss; the credible set has 60 SNPs at low purity. Explain why and recommend an HLA-specific workflow."

### Coloc Integration
> "Fine-map trait1 and trait2 separately at the same locus, then run coloc.susie. Report PP.H4 per credible-set pair."

### Reconciliation
> "SuSiE finds 3 credible sets but FINEMAP finds 1 at the same locus. Diagnose: is it convergence, LD mismatch, or non-sparse architecture?"

The workflow (locus window, LD source, LD diagnostic, method choice, credible-set reporting, coloc.susie hand-off) is the `SKILL.md` Decision Tree, Critical LD Diagnostic Block, Required Reporting Schema and Coloc.susie Integration.

Tips on LD reference choice, credible-set interpretation, L selection, purity filtering, PolyFun's `prior_weights` argument, non-sparse loci, HLA, PSD violations, and cross-ancestry gains are covered in `SKILL.md` (Per-Tool Failure Modes, Quantitative Thresholds, Cross-Ancestry Fine-Mapping with SuSiEx, and Common Errors).

## Related Skills

- causal-genomics/colocalization-analysis - coloc.susie operates on credible sets
- causal-genomics/mendelian-randomization - Fine-mapped cis-instruments
- causal-genomics/pleiotropy-detection - Per-credible-set pleiotropy
- population-genetics/linkage-disequilibrium - LD matrix construction
- population-genetics/association-testing - Upstream GWAS summary statistics
- workflows/gwas-pipeline - End-to-end GWAS to fine-mapping pipeline
- variant-calling/variant-annotation - Annotate credible-set variants
