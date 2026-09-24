# Proteome-Wide Drug-Target Mendelian Randomization - Usage Guide

## Overview

Cis-pQTL Mendelian randomization is the genetic-instrument analogue of a placebo-controlled drug trial: a genetic variant near the gene encoding a target protein perturbs that protein's plasma level from conception, and the variant's effect on a downstream phenotype estimates the causal effect of pharmacological inhibition of that target. The PCSK9-LDL-cholesterol-CAD story is the canonical positive example: cis-MR at the PCSK9 locus correctly predicted the direction and magnitude of clinical effect of evolocumab and alirocumab years before the FOURIER and ODYSSEY trials. The PCSK9 -> type 2 diabetes signal (Schmidt 2017 Lancet Diabetes Endocrinol) was the canonical on-target adverse-effect discovery.

This skill operationalises the Schmidt 2020 cis-MR framework (Nat Commun 11:3255) for plasma-proteome drug-target validation. It handles cis-pQTL instrument selection from UKB-PPP (Olink antibody), deCODE / Fenland / INTERVAL / ARIC (SomaScan aptamer), and FinnGen-PPP; cis-IVW with correlated instruments; colocalization triangulation; phenome-wide cis-MR for on-target adverse-effect scanning; cross-platform replication; and PAV (protein-altering-variant) sensitivity. The agent enforces the triangulation rule that a single significant cis-MR is necessary but never sufficient for a drug-target nomination.

## Prerequisites

R packages (TwoSampleMR, MendelianRandomization, coloc, susieR, ieugwasr, MR-PRESSO), Ensembl VEP, a local 1000 Genomes EUR plink reference and the pQTL summary-statistic sources are listed in SKILL.md under "Tool Installation Notes".

## Quick Start

Tell the AI agent what to do:
- "Run a cis-MR of PCSK9 on coronary artery disease using UKB-PPP cis-pQTLs"
- "Test IL6R inhibition on rheumatoid arthritis with cis-pQTL instruments from deCODE plus coloc triangulation"
- "Run a phenome-wide cis-MR of PCSK9 across OpenGWAS to find on-target adverse effects"
- "Replicate this cis-MR target across UKB-PPP Olink and deCODE SomaScan and flag PAV-confounded instruments"
- "Test a target with two independent cis-signals using coloc.susie per credible set"
- "Apply the Burgess 2016 sample-overlap correction to a UKB-PPP exposure with a UKB phenotype outcome"

## Example Prompts

### Single Drug Target, Single Outcome

> "I have UKB-PPP cis-pQTLs for PCSK9 within +/-500 kb of the gene and CARDIoGRAMplusC4D CAD summary stats. Run cis-IVW with weak-IV filtering, harmonise with action=2, triangulate with coloc.abf at p12=5e-6, and report PP.H4 plus PAV-excluded sensitivity. Annotate every cis-pQTL with VEP first."

> "Test IL6R protein -> rheumatoid arthritis using cis-pQTLs from deCODE SomaScan plus colocalization. Flag any cis-pQTL coloc'd with neighbouring genes' eQTLs in GTEx whole blood."

### Phenome-Wide On-Target Adverse-Effect Scan

> "Hold the PCSK9 cis-pQTL instrument set fixed and run cis-MR against all OpenGWAS outcomes with sample size >= 50,000 in European-ancestry cohorts. Bonferroni-correct over outcomes and report a phewas-style forest plot of significant hits."

> "Screen IL23R cis-pQTL effects across all FinnGen DF12 disease endpoints for on-target adverse-effect discovery prior to a clinical-trial protocol."

### Target with Multiple Independent Cis-Signals (Allelic Heterogeneity)

> "ANGPTL3 has two independent cis-pQTLs in low LD. Run coloc.susie on the cis-window with in-sample LD, report per-credible-set PP.H4 against triglycerides GWAS, and run a Wald ratio per credible set."

### Cross-Platform Replication

> "Run the same cis-MR of target X on disease Y independently in UKB-PPP (Olink) and deCODE (SomaScan). Report direction agreement, magnitude ratio, and flag the protein as platform-discordant if direction disagrees."

### Sample-Overlap Correction

> "Both my exposure GWAS (UKB-PPP cis-pQTL for protein X) and outcome GWAS (UKB HES-derived phenotype) are from UK Biobank. Apply MR-RAPS with the one-sample-equivalent treatment, or switch the outcome to FinnGen for an independent-cohort test."

### Drug Repurposing / Target Nomination

> "Cross-reference cis-MR estimates with Open Targets L2G scores and the Open Targets Drug platform to nominate druggable proteins for an autoimmune indication. Require cis-MR P < 1.7e-5, coloc PP.H4 >= 0.7, and cross-platform replication."

## What the Agent Will Do

It follows the Cis-MR Standard Workflow and the Triangulation Requirement in SKILL.md, and reports the claim-ladder rung the evidence reaches.

## Related Skills

causal-genomics/mendelian-randomization - Parent polygenic-MR framework; cis-MR is the drug-target specialization
causal-genomics/colocalization-analysis - Required PP.H4 triangulation for any cis-MR drug-target claim
causal-genomics/fine-mapping - Credible-set construction prior to coloc.susie at the cis-locus
causal-genomics/pleiotropy-detection - MR-PRESSO / Egger diagnostics adapted to cis-window
causal-genomics/transcriptome-wide-association - eQTL-based parallel evidence for the same target
causal-genomics/mediation-analysis - Step from cis-MR to downstream mediator pathway
population-genetics/association-testing - Source GWAS pipelines for pQTL discovery
population-genetics/linkage-disequilibrium - LD-matrix construction for correlated cis-IVW
variant-calling/variant-annotation - VEP PAV annotation for sensitivity analysis
clinical-databases/clinvar-lookup - Pathogenic-variant context for nominated targets
