# Effector Gene Prioritization Skill Usage Guide

## Overview

Effector gene prioritization is the central interpretation step that bridges variant-level statistical fine-mapping and gene-level biological hypothesis at GWAS loci. The agent takes a GWAS lead variant (or credible set), interrogates a portfolio of variant-to-gene (V2G) methods (Open Targets L2G, MAGMA, FUMA, cS2G, PoPS, FUSION/S-PrediXcan TWAS, ABC, ENCODE-rE2G), and integrates orthogonal evidence streams (fine-mapping PIP, colocalization PP.H4, distance, chromatin enhancer-gene linking, polygenic similarity priors) into a concordance-based confidence tier.

No single tool is sufficient for this call; see `SKILL.md`'s Decision Tree, Per-Method Failure Modes, and Multi-Evidence Integration Framework for the operational rules, thresholds, and pitfalls (nearest-gene assumption, tissue mis-specification, MAGMA window choice, PoPS/L2G discordance, multi-effector loci). This Skill is research-level and population-scale -- see `SKILL.md`'s Scope note before applying it to any individual-patient question.

## Prerequisites

See `SKILL.md`'s Prerequisites section for required inputs, the eQTL/pQTL panel selection table, and Tool Install Notes for install commands.

## Quick Start

Tell the AI agent what to do:
- "Prioritise effector genes at my GWAS lead locus using Open Targets L2G + PoPS + coloc concordance"
- "Run MAGMA gene-based and gene-set enrichment on my GWAS sumstats, then layer PoPS for polygenic priority"
- "I have a fine-mapped credible set; which gene is the most likely causal effector?"
- "Triangulate L2G, PoPS, coloc PP.H4, ABC, and distance to nominate a candidate causal gene at the SORT1 locus"
- "Run FUMA SNP2GENE on my GWAS and compare with Open Targets L2G calls for the top 20 hits"
- "Use ABC enhancer-gene predictions to assign a distal regulatory variant to its likely target gene"

## Example Prompts

### Open Targets covers the trait
> "My GWAS is on type 2 diabetes (UKB-2020). Query Open Targets Genetics for the L2G top-ranked genes at all genome-wide-significant loci. Cross-check each with V2G sub-scores and flag genes where yProbaDistance dominates yProbaMolecularQTL (distance-only candidates). Report the top 50 with L2G score >= 0.5."

### Custom trait from scratch
> "I have GWAS summary statistics for a custom rare-disease phenotype not in Open Targets. Run MAGMA gene-based with a 35kb upstream + 10kb downstream window using the 1000G EUR LD reference, then run PoPS on the MAGMA Z output with the pre-built feature matrix. Cross-reference with SuSiE fine-mapping credible sets and coloc PP.H4 against eQTLGen whole blood. Report concordance per gene."

### Tissue-known prioritisation
> "For my LDL-cholesterol GWAS, the causal tissue is liver. Run S-PrediXcan with GTEx liver MASHR weights and coloc.susie with GTEx liver eQTL panel. Layer ABC enhancer-gene predictions from HepG2 (ENCODE) for distal regulation. Integrate with Open Targets L2G and report effector-gene candidates with >= 3 concordant evidence streams."

### Tissue-unknown prioritisation
> "I have a schizophrenia GWAS but the causal cell type within brain is unclear. Run LDSC-SEG to prioritise brain tissues, then S-MultiXcan across all GTEx brain sub-regions, then PoPS for polygenic priority. Report per-gene concordance across L2G, PoPS, and tissue-prioritised coloc."

### Distal regulation suspected
> "At my chr8:9p21 GWAS lead, the nearest gene is CDKN2A but I suspect long-range regulation. Run ABC and ENCODE-rE2G in matched cell type (CMK or vascular smooth muscle), cross-check with Open Targets L2G, and report whether the fine-mapped credible variant maps to CDKN2A or a distal gene (e.g. MTAP, ANRIL)."

### Publication-grade triangulation
> "I am writing up an effector-gene nomination for ANGPTL4 at a triglycerides GWAS lead. Triangulate (a) fine-mapping (SuSiE PIP), (b) coloc PP.H4 with GTEx subcutaneous adipose eQTL, (c) Open Targets L2G score, (d) PoPS score, (e) ABC enhancer-gene linking in adipocytes (Engreitz portal), (f) distance to TSS. Report concordance per locus and flag whether ANGPTL4 meets the >= 3 of 6 high-confidence threshold."

### Multi-effector locus
> "At my chr11p15.5 lipid GWAS locus, three genes (CLU, NCAM1, MTHFD1L) all show modest L2G scores. Run conditional analysis (FUSION.post_process.R conditional/joint analysis, run by default, or GCTA-COJO) to test independence, then per-gene coloc.susie at each independent signal, and report whether the locus is multi-effector or LD-tied."

## Related Skills

- causal-genomics/fine-mapping - Variant-level credible sets feeding L2G and concordance scoring
- causal-genomics/colocalization-analysis - coloc PP.H4 as one of the six evidence streams
- causal-genomics/transcriptome-wide-association - Gene-level association and FOCUS gene fine-mapping
- causal-genomics/mendelian-randomization - cis-eQTL MR as orthogonal causal evidence
- causal-genomics/mediation-analysis - Downstream gene-mediated trait effects
- causal-genomics/proteome-mr-drug-target - pQTL-based drug-target nomination
- atac-seq/enhancer-gene-linking - ABC and ENCODE-rE2G enhancer-gene predictions
- atac-seq/atac-peak-calling - Enhancer candidates in matched tissue
- gene-regulatory-networks/coexpression-networks - Gene co-expression features feeding PoPS
- gene-regulatory-networks/scenic-regulons - TF-target regulons as supporting evidence
- pathway-analysis/go-enrichment - Pathway context for candidate effector genes
- population-genetics/association-testing - Upstream GWAS summary-statistic generation
- variant-calling/variant-annotation - Coding-consequence annotation
- workflows/gwas-pipeline - End-to-end GWAS pipeline producing input
