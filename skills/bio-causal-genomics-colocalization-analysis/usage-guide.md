# Colocalization Analysis - Usage Guide

## Overview

Colocalization tests whether two or more association signals at a genomic locus are driven by the same causal variant. It is the standard tool for integrating GWAS hits with molecular QTLs (eQTL, sQTL, pQTL, mQTL), prioritising candidate genes, and distinguishing shared causality from coincidental overlap due to linkage disequilibrium. Methods range from the fast single-causal coloc.abf (Giambartolomei 2014) to the multi-causal coloc.susie (Wallace 2021), the multi-trait HyPrColoc (Foley 2021), and the pleiotropy-vs-linkage SMR/HEIDI (Zhu 2016).

The agent will pick the appropriate method based on the experimental scenario: number of traits, expected number of causal variants per locus, availability of LD reference, ancestry composition, and whether the question is causal-mediation (SMR) or shared-variant (coloc). It will harmonise allele coding across summary statistics, format inputs correctly (`type`, `sdY`, `s`, `N`), run colocalization, and ALWAYS report sensitivity over the p12 prior alongside the headline PP.H4.

## Prerequisites

See SKILL.md "Tool Install Notes" for the exact install command per tool (coloc, susieR, HyPrColoc, SMR, eCAVIAR, PWCoCo, SharePro_coloc, moloc, plus the ggplot2/patchwork/data.table plotting deps). LD reference panel extraction via plink2 is in SKILL.md "LD Matrix Construction for coloc.susie".

## Quick Start

Tell the AI agent what is needed:
- "Test if my GWAS lead SNP at chr6:30500000 colocalizes with the IL2 eQTL in whole blood"
- "Run coloc.susie at this locus; I have two independent GWAS signals after conditional analysis"
- "Check colocalization between my GWAS and all 49 GTEx tissues using HyPrColoc"
- "I have a GWAS + eQTL + sQTL + mQTL at the same locus; run moloc"
- "Harmonise these summary stats and run coloc with sensitivity over p12"
- "Run SMR + HEIDI between my GWAS and eQTLGen blood eQTL"
- "The locus is in the MHC -- what method should be used"
- "My GWAS is in East Asian ancestry but the eQTL is GTEx EUR -- can coloc be trusted"

## Example Prompts

### Single-Tissue GWAS-eQTL

> "Take this 1 Mb window centred on the GWAS lead SNP and run coloc.abf against the eQTL for the nearest gene. Report PP.H4 with p12 sensitivity."

> "I have summary stats with p-values and MAF but no betas -- run coloc using p-value / MAF input format."

### Multi-Causal / Allelic Heterogeneity

> "GCTA-COJO conditional analysis identified two independent GWAS signals at this locus. Run coloc.susie with the LD matrix and report all credible-set pair PPs."

> "PP.H3 is dominating coloc.abf despite obvious overlap in LocusZoom. Switch to coloc.susie or eCAVIAR and re-test."

### Multi-Tissue

> "Run coloc.abf between this GWAS locus and the same gene across all 49 GTEx v8 tissues; rank tissues by PP.H4 to identify causal cell type."

> "Use HyPrColoc across GTEx tissues to cluster tissues that share the causal variant."

### Multi-Trait / Multi-Omic

> "Integrate GWAS + eQTL + sQTL + pQTL at this locus with moloc and report the PPA for the all-share hypothesis."

> "Cluster 12 cardiovascular GWAS traits at the LDLR locus with HyPrColoc."

### SMR / HEIDI

> "Run SMR + HEIDI between my disease GWAS .ma file and eQTLGen .besd. Report SMR p, HEIDI p, and number of HEIDI SNPs."

### MHC / HLA

> "This locus is in the MHC. Flag the long-range-LD problem and recommend HLA-coloc (Butler-Laporte 2024) or exclude-and-report-HLA-association strategy."

### Ancestry-Mismatched

> "GWAS is FinnGen (FIN) and eQTL is GTEx (EUR). Compute z-score vs LD diagnostics with estimate_s_rss; if lambda > 0.05, switch to coloc.abf or SharePro_coloc."

### Visualization

> "Make a LocusCompare plot and a stacked regional Manhattan for this colocalization."

## What the Agent Will Do

1. Identify the question type (single-causal, multi-causal, multi-trait, pleiotropy vs linkage) and select method from the decision tree.
2. Extract a +/- 500 kb to 1 Mb window centred on the GWAS lead SNP (or joint top-variant when both traits available).
3. Harmonise allele coding between summary statistics; flip betas where A1/A2 swap; drop palindromic SNPs at high MAF.
4. Format the per-trait coloc input lists (beta, varbeta, snp, position, N, type, sdY for quant, s for cc).
5. Run the chosen colocalization method (coloc.abf, coloc.susie, HyPrColoc, moloc, SMR, eCAVIAR, PWCoCo, or SharePro_coloc).
6. For coloc.susie: run `susieR::estimate_s_rss` to verify z-score vs LD consistency; abort if lambda > 0.05.
7. ALWAYS run `coloc::sensitivity(res, 'H4 > 0.75')` and report the p12 range over which PP.H4 stays above threshold.
8. Interpret PP.H0-H4 against the appropriate threshold (>= 0.75 screening, >= 0.80 published, >= 0.90 stringent).
9. Generate regional association and LocusCompare plots colored by LD to lead.
10. For multi-causal results, report per-(CS1, CS2) PP and lead SNP per credible set.
11. For multi-omic results, report the all-share PPA and per-pair PPs.
12. Flag known failure modes (MHC, ancestry mismatch, low-N eQTL, lead-SNP-swap window bias) explicitly in the report.

## Plain-Language H0-H4 (for Methods Section)

For a methods or supplementary description in plain prose:

- **H0** -- Neither trait has an association signal at this locus. Posterior reflects "no signal anywhere".
- **H1** -- Only trait 1 has a causal variant in the window; trait 2 has no signal.
- **H2** -- Only trait 2 has a causal variant in the window; trait 1 has no signal.
- **H3** -- Both traits have causal variants in the window, but they are different SNPs (linkage / coincidence under LD).
- **H4** -- Both traits share a single causal variant in the window (colocalization).

A high PP.H4 (>= 0.75) supports a shared-causal-variant interpretation. A high PP.H3 supports distinct causal variants in linkage. PP.H0 / PP.H1 / PP.H2 indicate the locus is underpowered for at least one trait. Always report all five posteriors, not PP.H4 alone.

Worked harmonisation code + pitfalls, the PWCoCo conditional CLI recipe, and the lead-SNP-swap operational steps now live in SKILL.md (Allele Harmonisation, PWCoCo, and Lead-SNP swap and window bias sections respectively) -- see there rather than duplicating here.

## Related Skills

- causal-genomics/mendelian-randomization - Downstream MR using colocalized SNPs as IVs
- causal-genomics/fine-mapping - SuSiE credible sets that feed coloc.susie
- causal-genomics/mediation-analysis - Causal mediation building on shared causal variants
- causal-genomics/pleiotropy-detection - Distinguishing horizontal pleiotropy from shared causality
- population-genetics/association-testing - GWAS summary stat generation
- population-genetics/linkage-disequilibrium - LD panel construction for coloc.susie / PWCoCo
- variant-calling/variant-annotation - Functional annotation for variant-specific priors
- single-cell/scatac-analysis - Per-cell-type chromatin context for coloc results
- differential-expression/deseq2-basics - eQTL count generation
- workflows/gwas-pipeline - Upstream GWAS producing coloc input
