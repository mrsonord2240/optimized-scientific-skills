# Library Design - Usage Guide

## Overview

Design pooled sgRNA libraries for CRISPR knockout, interference (CRISPRi), activation (CRISPRa), Cas12a multiplex, base-editor, or prime-editor screens. Covers chemistry selection, on-target scoring (Rule Set 2 / Azimuth / DeepSpCas9 / CRISPRon), off-target scoring (CFD / MIT / Elevation), TSS-relative positioning, control composition, oligo synthesis layout, and post-cloning QC.

## Quick Start

Tell the AI agent what to design:
- "Design a custom Brunello-style Cas9 knockout library for 250 kinases with 4 guides per gene plus 500 non-targeting controls"
- "Design a CRISPRi library against 1,200 lncRNAs using Dolcetto positioning rules against the FANTOM5 TSS"
- "Build a CRISPRa Calabrese-style library targeting -150 to -75 of TSS for 800 transcription factors"
- "Design an enAsCas12a paralog library covering 600 paralog pairs as 4-guide in4mer arrays"
- "Design a base-editor library tiling editing windows across all exons of BRCA1 and BRCA2 for variant scanning"
- "Diagnose why my freshly cloned plasmid pool has Gini 0.28 and ~3% zero-count guides"

## Example Prompts

### Cas9 Knockout Libraries

> "Design a focused Cas9 KO library targeting 350 DNA damage response genes. Use Rule Set 2 / Azimuth on-target scoring and CFD off-target filtering. Target the first 5-65% of each protein, exclude guides with GC outside 30-70% or poly-T runs, 4 guides per gene, 500 non-targeting controls, 50 AAVS1 controls, 50 CEGv2 reference essentials, 50 NEGv1 non-essentials."

> "I'm screening in HCT116 cells. Should I use DeepSpCas9 instead of Azimuth Rule Set 2 for guide ranking? Compare the two and recommend."

> "I have a custom gene list of 80 paralog pairs. Should I use Cas9 single-KO or Cas12a multiplex, and why?"

### CRISPRi / CRISPRa Libraries

> "Design a Dolcetto-style CRISPRi library for 600 transcription factors. Resolve TSS from FANTOM5 highest-rank CAGE peak per gene; if FANTOM5 lacks a peak, fall back to Ensembl canonical TSS and flag those genes for review. Search the Dolcetto window (-50 to +300 from TSS), ranking toward the +25 to +75 optimum. 6 guides per gene plus 1,000 NTCs."

> "Compare Calabrese vs Horlbeck CRISPRa TSS targeting windows for my activation screen in iPSC-derived neurons. We need to detect modest fold-change activation."

> "My CRISPRi screen has weak dropout signal even on RPL/RPS genes. Audit guide positioning against current FANTOM5 / matched neuron CAGE data and re-design any guides outside ±100 from the empirical TSS."

### Cas12a Multiplex / Paralog Libraries

> "Build an Inzolia-style enAsCas12a 4-guide array library covering 400 paralog pairs in the receptor tyrosine kinase family. Include singleton controls (gene A alone, gene B alone, double-NTC) so we can score genetic interaction = double_KO_LFC - sum(single_KO_LFC)."

> "Design an in4mer triple-KO library testing all triplets within a 25-gene synthetic-lethality hypothesis."

### Base / Prime Editor Libraries

> "Design a CBE saturation library across BRCA1 exons 1-23 tiling every NGG-adjacent spacer that places at least one C in editing positions 4-8. Flag bystander Cs and annotate predicted amino acid change."

> "Build a prime-editor library to install 320 specific ClinVar variants in MLH1, MSH2, MSH6, PMS2. Use PRIDICT2 to score pegRNA candidates and select the top-3 scored pegRNAs per intended edit."

### Library Diagnostics

> "My plasmid pool shows Gini 0.28, skew ratio 6, and 2.5% zero-count guides on Brunello-style sequencing. Diagnose the cause (PCR bias, synthesis defect, cloning bottleneck, library age) and recommend remediation."

> "We see ERBB2 dropping out as 'essential' in HER2-amplified SK-BR-3 cells. Is this a library-design issue or a copy-number bias problem?"

## What the Agent Will Do

The design workflow, thresholds, controls, oligo layout and QC targets are in `SKILL.md` (start at "sgRNA Library Design": it lists the required inputs the agent confirms before designing, and the deliverable).

## Related Skills

- crispr-screens/screen-qc - Library skew, Gini, replicate correlation, essentialome PR-AUC
- crispr-screens/mageck-analysis - Standard analysis pipeline for the designed library
- crispr-screens/combinatorial-screens - Cas12a multiplex / paralog-pair library design
- crispr-screens/base-editing-analysis - base-editor library design and editing-window analysis
- crispr-screens/prime-editing-screens - PRIDICT2-optimized pegRNA libraries
- crispr-screens/copy-number-correction - Filter amplicon-driven artifacts in cancer-cell-line screens
- crispr-screens/in-vivo-screens - Bottleneck math for animal-model focused libraries
