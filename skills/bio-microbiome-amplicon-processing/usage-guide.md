# Amplicon Processing - Usage Guide

## Overview

This skill turns demultiplexed marker-gene FASTQ (16S rRNA, ITS) into an ASV (amplicon sequence variant) feature table plus representative sequences, using DADA2's per-run error model. An ASV is a model-inferred exact sequence conditioned on one sequencing run - not a clustered OTU and not an organism. The decisions that matter (which are invisible in a recipe) are: remove primers before truncating, learn the error model separately per run, choose truncation lengths within the merge-overlap budget, and never fix-truncate variable-length ITS. The output feeds taxonomy assignment, diversity, and differential abundance; the compositional statistics of the resulting table are shared with shotgun metagenomics.

Install commands, the method, thresholds and failure modes are in `SKILL.md`. You need demultiplexed FASTQ, the primer sequences and the amplicon region; for low-biomass samples, also negative controls.

## Quick Start

Tell your AI agent what you want to do:
- "Remove primers from my paired-end 16S reads, then run DADA2 to get an ASV table"
- "I have V3-V4 reads on 2x250 and my merge rate is near zero - help"
- "Process two MiSeq runs correctly with a per-run error model and merge the tables"
- "Process my fungal ITS amplicons into ASVs"
- "These are low-biomass skin swabs with extraction blanks - remove reagent contaminants with decontam"

## Example Prompts

### ASV inference
> "I have demultiplexed paired-end 16S V4 reads. Remove the primers with cutadapt, learn the error model per run, denoise, merge pairs, and remove chimeras to give me an ASV table and a read-tracking summary."

### Truncation budget
> "My amplicon is V3-V4 (~460 bp) on 2x250 reads. What truncLen should I use, and how do I keep enough overlap to merge while controlling quality?"

### Multi-run studies
> "My samples came off three sequencing runs. Show me how to learn the error model on each run separately and combine the sequence tables before chimera removal."

### ITS / non-16S
> "I have fungal ITS2 amplicons of variable length. Trim the spacer and infer ASVs without fixed truncation."

### Platform caveats
> "These are NovaSeq reads with binned quality scores. Check whether the DADA2 error model is fit correctly and fix it if not."

### Low-biomass / decontamination
> "These are low-biomass biopsy samples sequenced with extraction-blank and no-template-PCR negative controls. After building the ASV table, run decontam to flag and remove reagent/kit contaminants, and report how many ASVs and reads were removed."

## Related Skills

- taxonomy-assignment - Assign taxonomy to the ASVs produced here
- diversity-analysis - Alpha/beta diversity of the resulting community table
- differential-abundance - Compositional DA on the ASV/feature table
- qiime2-workflow - The QIIME2 CLI equivalent of this R workflow
- read-qc/adapter-trimming - cutadapt primer removal before DADA2
- metagenomics/kraken-classification - Shotgun (not amplicon) read classification
- metagenomics/abundance-estimation - Shared compositional/normalization theory
- phylogenetics/tree-io - Phylogenetic tree for UniFrac / Faith PD
- workflows/microbiome-pipeline - End-to-end amplicon pipeline
