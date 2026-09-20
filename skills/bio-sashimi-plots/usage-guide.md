# Sashimi Plots - Usage Guide

## Overview
Visualize splicing events as sashimi-style plots showing per-sample read coverage and splice junction arcs labeled by read counts. Tools differ in input handling, group aggregation, and customization: ggsashimi (publication overlays), rmats2sashimiplot (rMATS-aware), MAJIQ-VOILA (LSV posteriors, licence-gated), leafviz (leafcutter Shiny), Jutils (tool-agnostic), pyGenomeTracks (multi-track figures). Install notes, flags, recipes and failure modes are in `SKILL.md`.

## Quick Start
Tell your AI agent what you want to do:
- "Create sashimi plots for the top 20 differential splicing events from rMATS"
- "Visualize a specific exon-skipping event with samples grouped by condition"
- "Generate publication-quality sashimi with intron shrinking and per-condition aggregation"
- "Plot splicing alongside ChIP-seq tracks for the same locus"

## Example Prompts

### Single Event Visualization
> "Plot a sashimi for chr17:43094000-43125000 (BRCA1 region) with control vs treatment groups, intron shrinking, and matched y-axis scales."

### Batch Plotting from Differential Output
> "Iterate ggsashimi over the top 25 significant rMATS events; output per-event PDFs with flanking 500nt context."

### MAJIQ VOILA
> "Run voila on MAJIQ deltapsi output to browse LSV posteriors and splice graphs." (needs a licensed MAJIQ install)

### leafviz
> "Generate leafcutter Shiny app for interactive cluster-level browsing."

### Multi-Track Figures
> "Use pyGenomeTracks to combine RNA-seq coverage tracks with H3K4me3 ChIP-seq for the same locus."

### Tool-Agnostic Heatmaps
> "Use Jutils to create a unified heatmap of significant events across rMATS and leafcutter output."

## Related Skills

- differential-splicing - Identify events to plot
- splicing-quantification - Context for PSI values
- data-visualization/genome-tracks - Multi-track figure design
- data-visualization/ggplot2-fundamentals - ggsashimi customization
