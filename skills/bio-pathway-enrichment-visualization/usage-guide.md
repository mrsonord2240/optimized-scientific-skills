# Enrichment Visualization - Usage Guide

## Overview
This skill turns an enrichment result object (an enrichResult from ORA, a gseaResult from GSEA, or a compareClusterResult) into a figure with the enrichplot package. The central decision is not which plotting function to call but how to handle gene-set REDUNDANCY: a default top-20 GO dotplot is usually one biological theme drawn twenty times, because the GO DAG and nested pathway databases guarantee that a real signal surfaces as a cluster of near-identical overlapping terms. The figure either SHOWS that redundancy as structure (emapplot/treeplot via pairwise_termsim, or EnrichmentMap) or DELETES it (simplify/REVIGO). The skill also covers keeping the NES sign for GSEA, distinguishing GeneRatio from fold enrichment, and admitting that showCategory truncates.

Install commands, version notes, and failure modes (including the `ggridges`/`ggarchery`/`ggupset` Suggests-only dependencies and the `treeplot`-on-`compareClusterResult` crash) live in SKILL.md's Prerequisites, Version Compatibility, and Common Errors sections - read those before running anything.

## Quick Start
Tell your AI agent what you want to do:
- "Make a dotplot of my GO enrichment, collapsing redundant terms first"
- "Show the redundant terms as clusters with an enrichment map"
- "Plot my GSEA results keeping the direction (activated vs suppressed)"
- "Show a GSEA running-score plot for my top pathway"

## Example Prompts

### Collapsing redundancy
> "My enrichGO BP result has 180 significant terms and the top-20 dotplot is full of cell-cycle synonyms. Collapse the GO-DAG redundancy with simplify(), then dotplot the top 20, and tell me how many terms survived."

> "I have an enrichResult with many overlapping terms. Run pairwise_termsim and draw an enrichment map so I can see which terms are really one biological theme, and a treeplot with five named clusters."

### Encoding and effect size
> "Make a dotplot of my GO results ordered by fold enrichment instead of GeneRatio, and explain why the two orderings differ."

> "Build a gene-concept network for my top 6 enriched terms colored by log2 fold change so I can see which genes bridge multiple terms."

### GSEA plots
> "Plot my gseaResult as a ridgeplot showing the leading-edge distribution per set, keeping direction, and a gseaplot2 for the single most significant pathway."

> "Show the running enrichment score for my top three GSEA pathways overlaid in one panel with the stats table."

### Comparison and saving
> "I ran compareCluster across up- and down-regulated genes. Make a faceted dotplot comparing the two."

> "Save my enrichment dotplot as a publication-quality PDF with a viridis color scale and a caption noting the total significant term count."

The agent's workflow (object-class identification, the redundancy decision, encoding choices, failure modes) and the full plot/class/redundancy reference table live in SKILL.md - see its Tool Taxonomy and Decision Tree sections.

## Related Skills

- go-enrichment - Produces the enrichResult; owns simplify()
- gsea - Produces the gseaResult; owns the enrichment score and leading edge
- kegg-pathways - KEGG results to plot
- reactome-pathways - Reactome results to plot
- wikipathways - WikiPathways results to plot
- data-visualization/ggplot2-fundamentals - Generic ggplot2 grammar for the returned objects
- workflows/expression-to-pathways - End-to-end DE-to-enrichment-to-figure pipeline
