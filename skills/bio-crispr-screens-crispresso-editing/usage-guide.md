# CRISPResso2 Editing Analysis - Usage Guide

## Overview

Decision-grade quantification of CRISPR editing outcomes with CRISPResso2 across Cas9 nuclease (indels + HDR), cytosine and adenine base editors (target conversion + bystander), and prime editor (templated edits) modes. Covers single-amplicon (`CRISPResso`), batch (`CRISPRessoBatch`), pooled (`CRISPRessoPooled`), WGS off-target (`CRISPRessoWGS`), and comparison (`CRISPRessoCompare`) workflows. Defines the quantification window math, the substitution-vs-indel ratio for distinguishing BE from Cas9 contamination, and the MMEJ deletion signature in allele tables.

## Prerequisites

Install and required inputs are in SKILL.md (Version Compatibility).

## Quick Start

Tell the AI agent what to analyze:
- "Quantify indel rate at my Cas9 cut site from amplicon sequencing"
- "Compute target editing % and bystander rate from my CBE experiment using quantification window size 10"
- "Run CRISPRessoBatch on 24 samples (timecourse: 0/6/12/24/48 hours) targeting BRCA1 exon 11"
- "Diagnose why my CRISPResso run has 35% alignment rate"
- "Distinguish BE-mediated edits from Cas9-contamination indels using substitution-vs-indel ratio"
- "Run CRISPRessoPooled on 96 arrayed-validation amplicons in a single MiSeq library"

## Example Prompts

### Cas9 Editing

> "Run CRISPResso on a single Cas9 sample at MLH1 exon 12. Amplicon: ACGT... (primer-trimmed). Guide: GUIDE_SEQ. Report % unmodified vs NHEJ; flag if indel rate <70% (incomplete KO)."

> "I have 50 samples from a Cas9 timecourse experiment in K562. Use CRISPRessoBatch with batch_settings.txt; parallelize across 16 processes. Output the aggregated CRISPRessoBatch_quantification_of_editing_frequency.txt for downstream plotting."

### Base Editor Analysis

> "Analyze my CBE sample at BRCA1 c.5135C>T target. CRISPResso with `--base_editor_output`, `--conversion_nuc_from C --conversion_nuc_to T`, quantification window size 10, center -10. Report: target conversion %, bystander rate at adjacent Cs, indel byproduct %, substitution-vs-indel ratio."

> "Compare ABE7.10 vs ABE8.20 editing efficiency at the same target site across 5 cell lines. Use CRISPRessoBatch; output a heatmap of % A->G at the target position."

> "Distinguish my BE3 sample's edits from background Cas9 indels: confirm substitution-vs-indel ratio >10 (clean BE) or <3 (Cas9-like cut activity); recommend whether to repeat with a tighter nCas9-BE3 vector."

### Prime Editor Analysis

> "Quantify PE-2 editing for installation of MLH1 c.677A>G. Specify spacer + extension (RTT+PBS) + scaffold sequences. Report intended-edit %, scaffold incorporation %, indel %; flag if scaffold incorporation >5% (RTT design issue)."

> "Compare 8 candidate pegRNAs for the same intended edit. Use CRISPRessoBatch to run on saturating-dose samples; report which pegRNA achieves highest intended/(scaffold+indel) ratio."

### Pooled-Amplicon Mode

> "Run CRISPRessoPooled on a 96-amplicon arrayed-validation pool. Verify amplicon-misassignment is <2% by checking the primer-overlap diagnostic."

### Comparisons

> "Compare CRISPResso outputs for control vs DNA-damage-treated samples using CRISPRessoCompare. Quantify shift in indel size distribution and MMEJ-like deletion enrichment."

### Diagnostics

> "My CRISPResso alignment rate is 38%. Diagnose: wrong amplicon sequence, primer-dimer contamination, low-quality FASTQ, or strand orientation issue?"

> "Allele table shows a -7 bp deletion in 35% of reads at the cut site. Distinguish MMEJ-mediated repair from random NHEJ; recommend follow-up validation."

## Related Skills

- crispr-screens/base-editing-analysis - Variant-function analysis using CRISPResso2 BE outputs
- crispr-screens/prime-editing-screens - PRIDICT2 pegRNA design + PE-tiling
- crispr-screens/library-design - sgRNA / pegRNA design rules
- crispr-screens/screen-qc - Editing-efficiency QC threshold for variant interpretation
- variant-calling/variant-annotation - Annotate edited variants downstream
- read-alignment/bwa-alignment - For WGS off-target alignment input
