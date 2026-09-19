# Taxonomy Assignment - Usage Guide

## Overview

Taxonomy assignment classifies amplicon ASVs/OTUs (16S, ITS, 18S) against a reference database. The output label is a classifier + database + primer-region-conditioned hypothesis at a stated confidence, NOT an identification. The single most important honest message: a short 16S read (one hypervariable region, ~250 bp for V4) licenses GENUS at best - frequently only family for poorly resolved clades - and rarely species. The confidence/bootstrap score measures the model's certainty assuming the true taxon is in the database; it cannot detect when the true taxon is absent, in which case the classifier returns the nearest wrong relative. ITS (the fungal barcode) is the exception and does resolve to species; 18S resolves coarsely like 16S.

Three classifier families are covered: DADA2 `assignTaxonomy` + `addSpecies` (RDP naive Bayes, R), DECIPHER IDTAXA (conservative tree-descent, R), and QIIME2 q2-feature-classifier (`classify-sklearn` naive Bayes and `classify-consensus-vsearch` alignment-consensus, CLI). Shotgun read classification is a different category - see metagenomics/kraken-classification and metagenomics/metaphlan-profiling.

## Prerequisites

```bash
# R classifiers (DADA2 + DECIPHER IDTAXA)
# BiocManager::install(c('dada2', 'DECIPHER'))

# QIIME2 installs as its own conda environment; the release tag defines the plugin API and the
# .qza classifier format (e.g. qiime2-amplicon-2024.10). vsearch ships with the QIIME2 env.
# conda env create -n qiime2-amplicon-2024.10 --file <release env file>
```

Conceptual prerequisites: inputs are per-feature sequences (a DADA2 `seqtab_nochim` or a QIIME2
`FeatureData[Sequence]`) - ASV inference happens upstream in amplicon-processing. See SKILL.md's
"Single Most Important Modern Insight" and Version Compatibility sections for what governs a
correct classification (database/region matching, release pinning).

## Quick Start

Tell your AI agent what you want to do:
- "Assign taxonomy to my 16S V4 ASVs against SILVA and report the rank the data supports"
- "Train a region-matched naive-Bayes classifier for my 515F/806R primers, then classify"
- "Classify my fungal ITS sequences against UNITE"
- "I am getting a scikit-learn version error on a pre-trained classifier - what are my options?"
- "Use IDTAXA for a conservative, novelty-aware classification"
- "Filter host mitochondria and chloroplast features out of my table before diversity and DA"

## Example Prompts

### Region-matched classification
> "I have 16S V4 (515F/806R) ASVs from DADA2. Assign taxonomy with SILVA, but make sure the classifier is matched to the V4 region rather than full-length, and report genus where species is not supported."

### Database selection
> "These are environmental samples with many under-named bacteria. Should I use SILVA or GTDB, and what changes about the names if I switch?"

### Method selection
> "Classify my ASVs two ways - naive Bayes and vsearch alignment-consensus - and tell me where they disagree and why one is immune to the scikit-learn version error."

### Confidence and honest reporting
> "Run assignTaxonomy and tell me what minBoot you used, what it trades, and which ASVs you left unassigned at genus rather than force-filling."

### Species-level question
> "Can I report species for this 16S ASV? If not, what would it take to defend a species call?"

### Fungal / eukaryote markers
> "My amplicon is fungal ITS. Classify against UNITE and explain why ITS can reach species when my 16S could not."

### Filtering host organelle reads
> "My samples are plant-associated and a big fraction of ASVs are labelled Chloroplast or Mitochondria. Filter the host organelle features out of the feature table after assignment, before diversity and differential abundance, and tell me what fraction of reads that removed."

## Related Skills

- amplicon-processing - Generate the ASV table that is classified here
- diversity-analysis - Alpha/beta diversity of the classified community table
- differential-abundance - Compositional differential abundance on the classified feature table
- qiime2-workflow - The QIIME2 CLI workflow this classification step plugs into
- read-qc/adapter-trimming - cutadapt primer removal before ASV inference and assignment
- metagenomics/kraken-classification - Shotgun (raw-read, not ASV) k-mer classification
- metagenomics/metaphlan-profiling - Shotgun marker-gene profiling; a different input artifact
- phylogenetics/tree-io - Phylogenetic tree for UniFrac / Faith PD on the classified table
- workflows/microbiome-pipeline - End-to-end amplicon pipeline
