# Transcriptome-Wide Association (TWAS) Skill Usage Guide

## Overview

Transcriptome-wide association studies (TWAS) test gene-level association with a GWAS trait via genetically predicted tissue expression. The agent runs TWAS from GWAS summary statistics using pre-trained eQTL prediction models (PrediXcan / FUSION / MetaXcan family), aggregates multi-tissue signals (S-MultiXcan / UTMOST), and probabilistically fine-maps co-significant gene clusters (FOCUS / MA-FOCUS).

TWAS is the natural cousin of Mendelian randomization and colocalization; a single TWAS hit is associational, not causal. SKILL.md holds the failure modes and the triangulation rule.

## Prerequisites

- GWAS summary statistics with standard columns (SNP, A1/A2, BETA or Z, SE, P, N or use_n)
- Pre-trained prediction weights: FUSION panels from gusevlab.org/projects/fusion (.pos + per-gene RData) OR PredictDB models from predictdb.org (GTEx v8 elastic-net or MASHR, .db + covariance)
- Ancestry-matched LD reference panel (e.g. 1000 Genomes EUR; PLINK .bim/.bed/.fam by chromosome)
- Python 3.9-3.11 + R 4.3+ + PLINK 1.9 + PLINK 2.0
- Disk space: GTEx v8 PredictDB MASHR-EUR is approximately 8 GB; FUSION GTEx weights per tissue are 1-3 GB

Install commands, pins and the required pyfocus patch are in SKILL.md, Tool Install Notes.

## Example Prompts

### Single-tissue TWAS from GWAS summary statistics
> "Run S-PrediXcan on my LDL-cholesterol GWAS using the GTEx v8 MASHR liver model. Filter to genes with significant association at p < 2.3e-6 (Bonferroni for 22k genes) and report the top 20 by Z magnitude. Note the assumption that EUR-trained weights apply to an EUR-only GWAS."

### Multi-tissue prioritisation
> "Run S-PrediXcan across all 49 GTEx v8 MASHR tissues for my schizophrenia GWAS. Combine the per-tissue outputs with S-MultiXcan for a joint multi-tissue test. Compare the joint-significant gene set against the per-tissue Bonferroni-significant lists, and flag any genes where the joint test gains power over the best per-tissue test."

### TWAS fine-mapping at a gene-dense locus
> "At the chr11p15.5 locus my TWAS reports 7 genes passing significance. Run FOCUS using the GTEx v8 whole blood DB and 1000G EUR LD reference to compute per-gene PIPs and report the credible gene set. Genes with PIP >= 0.8 are causal candidates; co-significant genes with PIP < 0.5 are LD-tagged."

### Multi-ancestry TWAS
> "Run MA-FOCUS on my T2D GWAS combining EUR (BBJ + UKB), EAS (BBJ), and AFR (AAGILE) sumstats using ancestry-matched 1000 Genomes LD references. Use the MA-FOCUS DBs for each ancestry's PredictDB v8 MASHR adipose tissue. Report joint PIPs and per-ancestry contribution."

### Triangulation with MR and coloc
> "I have a TWAS hit for SORT1 in liver from S-PrediXcan. Triangulate this with (a) a cis-eQTL MR using SORT1 cis-eQTLs from GTEx liver as the exposure and my LDL GWAS as the outcome with TwoSampleMR, and (b) coloc.abf between the GTEx liver cis-eQTL for SORT1 and the GWAS at the SORT1 locus. Report a 4-way concordance summary across TWAS, FOCUS PIP, cis-MR, and coloc PP.H4."

### Splice-TWAS for a neuropsychiatric trait
> "Run TWAS using GTEx splicing models (sQTL-based PredictDB) instead of expression for my major depressive disorder GWAS in brain frontal cortex. Identify splice-mediated gene hits and contrast with the expression-TWAS hit list at the same loci."

### FUSION conditional analysis
> "Run FUSION on my CAD GWAS with GTEx artery coronary weights for all 22 autosomes. Then run FUSION.post_process.R at every significant locus to identify conditionally independent genes; report the joint-Z table per locus."

## Related Skills

- causal-genomics/fine-mapping - Variant-level credible sets are the precursor to FOCUS gene-level fine-mapping
- causal-genomics/colocalization-analysis - Coloc PP.H4 triangulation with TWAS hits
- causal-genomics/mendelian-randomization - cis-eQTL MR triangulation; drug-target prioritisation
- causal-genomics/mediation-analysis - Downstream gene-mediated trait effects given TWAS hits
- population-genetics/association-testing - Upstream GWAS summary statistic generation
- population-genetics/linkage-disequilibrium - LD reference panel construction
- differential-expression/deseq2-basics - eQTL count data for custom prediction-weight training
- single-cell/preprocessing - Cell-type-resolved eQTL panels for sc-TWAS
- workflows/gwas-pipeline - Upstream GWAS pipeline producing TWAS input
- variant-calling/variant-annotation - Functional annotation of TWAS / FOCUS top variants
