# Diversity Analysis - Usage Guide

## Overview

Diversity analysis summarizes the whole microbial community: alpha diversity (within-sample richness and evenness) and beta diversity (between-sample dissimilarity) of an amplicon ASV/OTU table, with QIIME2 core-metrics-phylogenetic, phyloseq/vegan, and scikit-bio. The central discipline is that a diversity number is the output of three choices made before the number appears - the rarefaction depth, the tree, and the metric - so the agent declares all three and shows the conclusion survives a second reasonable choice rather than reading a single p-value.

The two decisions that dominate everything: the sampling depth in `core-metrics` silently deletes every sample below it (biased toward low-biomass samples), and the tree behind UniFrac/Faith PD is a modeling choice (de novo from short reads loses to SEPP fragment-insertion and Greengenes2 placement). Per-taxon testing belongs to differential-abundance; shotgun profiler tables belong to metagenomics/metagenome-visualization; the shared compositional and rarefaction-debate theory lives in metagenomics/abundance-estimation, and the Hill-number and PERMANOVA-dispersion theory in metagenomics/metagenome-visualization.

## Prerequisites

See SKILL.md's "Required Setup" for install commands and data prerequisites.

## Quick Start

Tell your AI agent what you want to do:
- "Pick a rarefaction sampling depth from my feature table and tell me which samples it drops"
- "Run core-metrics-phylogenetic and report alpha and beta diversity"
- "Calculate Shannon and Faith PD per sample and compare my groups"
- "Compute weighted and unweighted UniFrac and test the group difference with PERMANOVA and betadisper"
- "Build a SEPP fragment-insertion tree instead of a de novo tree for UniFrac"

## Example Prompts

### Choosing the sampling depth
> "I have a QIIME2 feature table. Summarize the per-sample frequencies, show me an alpha-rarefaction curve, recommend a sampling depth on the plateau, and tell me how many and which samples that depth would drop."

### Alpha diversity
> "Rarefy my ASV table to the depth I chose, calculate observed features, Shannon, and inverse Simpson, report effective species, and test whether diversity differs between control and treatment with a non-parametric test."

### The tree decision
> "My reads are V4 16S. Build a SEPP fragment-insertion tree against a full-length reference instead of a de novo MAFFT+FastTree tree, filter out the ASVs that failed to insert, and explain why this matters for UniFrac."

### Beta diversity
> "Compute weighted and unweighted UniFrac plus generalized UniFrac at alpha 0.5, make PCoA plots, run PERMANOVA on each, and pair every PERMANOVA with betadisper so I know whether a significant result is a location shift or a dispersion difference."

### Compositional ordination
> "Run a robust Aitchison PCA (RPCA) on my ASV table and tell me which ASVs load on the first axis."

## Related Skills

- amplicon-processing - Generate the ASV table and representative sequences upstream
- taxonomy-assignment - Label the ASVs summarized here
- differential-abundance - Per-taxon between-group testing on the unrarefied counts
- qiime2-workflow - The QIIME2 CLI home for core-metrics and tree building
- phylogenetics/tree-io - Read, write, and root the UniFrac/Faith PD tree
- metagenomics/abundance-estimation - Shared CoDA and rarefaction-debate theory
- metagenomics/metagenome-visualization - Shared Hill-number and PERMANOVA-dispersion theory; diversity/ordination on shotgun profiler tables
- data-visualization/ggplot2-fundamentals - Custom ordination and diversity plots
