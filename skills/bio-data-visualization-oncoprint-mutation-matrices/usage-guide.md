# OncoPrint and Mutation Matrix Plots - Usage Guide

## Overview

Use this Skill to turn somatic-variant cohorts into multi-class gene-by-sample
OncoPrints, add clinical or burden tracks, and inspect co-occurring or mutually
exclusive driver pairs. It supports ComplexHeatmap for cohort-complete custom
plots, maftools for a rapid MAF view, and comut.py for Python workflows.

The `SKILL.md` sections “Build a cohort-complete matrix first” and “Common
Errors” contain the required class map, denominator rules, ID-alignment guards,
ordering behavior, and tool-specific caveats.

## Example Prompts

- “Build a cohort-complete OncoPrint of the top 20 recurrently mutated genes,
  stacked by alteration class, with subtype and stage annotations.”
- “Split the OncoPrint columns by molecular subtype and state that the displayed
  gene percentages remain cohort-wide.”
- “Order genes by the number of mutated samples and samples by TMB.”
- “Run somaticInteractions on the top 20 genes and return its pairwise result
  table with raw 2-by-2 counts, adjusted p-values, and effect directions.”
- “Create the equivalent burden-sorted comut.py plot while retaining cohort
  members with no mutation in the displayed genes.”
- “Log-transform the TMB track so hypermutators do not dominate.”

## Related Skills

- data-visualization/heatmaps-clustering - generic heatmaps
- data-visualization/lollipop-protein-maps - per-gene protein maps
- data-visualization/color-palettes - alteration-class palettes
- clinical-databases/variant-prioritization - upstream filtering
- variant-calling/variant-annotation - consequence annotation
- copy-number/cnv-annotation - CNV calls for OncoPrint cells
