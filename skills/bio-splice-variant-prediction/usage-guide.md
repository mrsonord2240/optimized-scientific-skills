# Splice Variant Prediction - Usage Guide

## Overview
Predict whether a DNA variant alters mRNA splicing with SpliceAI, Pangolin, MMSplice, SpliceTransformer, CI-SpliceAI, SpliceVault and CADD-Splice, and label the result with the ClinGen SVI 2023 PP3/BP4 evidence codes. Research-use decision support for expert review, not a diagnosis. Installation, commands, thresholds, failure modes and the Pangolin masking trade-off are in `SKILL.md` (sections Install, SpliceAI Workflow, Pangolin, Concordance).

## Quick Start
Tell your AI agent what you want to do:
- "Predict SpliceAI delta scores for variants in my VCF"
- "Score a panel of variants for splice impact with concordance across SpliceAI, Pangolin, and MMSplice"
- "Identify deep-intronic pseudoexon-creating variants using extended-window SpliceAI"
- "Label a candidate variant with ClinGen SVI 2023 splicing evidence (PVS1/PP3/BP4) for expert review"
- "Predict tissue-specific splicing impact with Pangolin for a brain-disease variant"

## Example Prompts

### Single Variant
> "For TP53 c.673-2A>G (GRCh38 chr17:7674292 T>C), compute SpliceAI delta scores with -D 50 and the Pangolin per-tissue scores (heart, liver, brain, testis). Check REF against my FASTA first."

### VCF Annotation
> "Run SpliceAI on my VCF with the grch38 reference, then label PP3/BP4 with the ClinGen thresholds and list any input records SpliceAI did not score."

### Deep-Intronic
> "I have an unsolved Mendelian case; re-run SpliceAI with -D 500 on my intronic candidates and report both pseudoexon boundaries."

### Concordance
> "Run SpliceAI + Pangolin + MMSplice on a candidate VUS panel; flag discordant predictions and variants a tool skipped."

### ASO Design
> "Walk me through the design considerations for a splice-switching ASO targeting an ESE region (checklist only)."

### Empirical Outcome
> "For a canonical 5'ss-disrupting variant, query SpliceVault to predict whether it causes exon skipping vs cryptic site activation."

## Related Skills

- splicing-qc - MaxEntScan + library QC for confirming predicted impact
- splicing-quantification - Empirical PSI from RNA-seq to validate predictions
- outlier-splicing-detection - FRASER2/DROP for RNA-seq confirmation in clinical samples
- variant-calling/clinical-interpretation - Broader ACMG/AMP variant interpretation
- variant-calling/variant-annotation - VEP plugin integration for SpliceAI
