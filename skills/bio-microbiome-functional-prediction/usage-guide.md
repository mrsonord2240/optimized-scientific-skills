# Functional Prediction - Usage Guide

## Overview

PICRUSt2 predicts community functional POTENTIAL from 16S amplicon ASVs by phylogenetic interpolation: it places each ASV on a tree of ~20,000 reference genomes and reports the gene content of its nearest sequenced relatives as the community's KO/EC/MetaCyc profile. The single fact that governs every decision: this is PREDICTED potential inferred from who-is-there, never measured gene content and never activity or expression - see SKILL.md's "Single Most Important Modern Insight" and "Predicted < Measured Ladder" sections. Accuracy is entirely a function of how well the community is represented in the reference set - credible in the densely-referenced human gut, collapsing toward a restatement of taxonomy in soil, marine, and novel environments.

Install commands and the full pipeline invocation are in SKILL.md's "Run the Pipeline" section.

## Quick Start

Tell your AI agent what you want to do:
- "Predict KO and MetaCyc potential from my 16S ASV table with PICRUSt2"
- "Run PICRUSt2 and report the NSTI distribution and the fraction of reads dropped"
- "Tell me whether PICRUSt2 is appropriate for my soil samples or whether I should use FAPROTAX"

## Example Prompts

### Prediction
> "I have a representative-sequence FASTA and an ASV abundance table from DADA2. Run PICRUSt2 with the recommended maximum-parsimony hidden-state method and produce KO, EC, and MetaCyc pathway tables."

> "Run the q2-picrust2 full pipeline on my QIIME2 feature table and rep-seqs."

### Quality gating
> "Summarize the NSTI distribution from my PICRUSt2 run and report how many ASVs and what fraction of reads were dropped at the default NSTI cutoff."

> "My mean NSTI is high - is this prediction trustworthy for marine sediment, and what should I use instead?"

### Method choice
> "I want to know whether this community is nitrifying. Should I use PICRUSt2 or FAPROTAX?"

> "I need measured functional gene content, not a prediction - what should I do?"

### Downstream
> "Run compositionally-aware differential abundance on my predicted MetaCyc pathways across two groups and frame the result correctly."

## Related Skills

- amplicon-processing - Generate the ASV table and representative sequences consumed here
- taxonomy-assignment - Taxonomic labels for the same ASVs
- differential-abundance - Compositional DA of the predicted KO/pathway table
- qiime2-workflow - The q2-picrust2 plugin path inside QIIME2
- metagenomics/functional-profiling - MEASURED shotgun function (HUMAnN); the predicted-vs-measured wall
- pathway-analysis/go-enrichment - Reading/enriching the predicted KO/MetaCyc lists
- workflows/microbiome-pipeline - End-to-end amplicon pipeline
