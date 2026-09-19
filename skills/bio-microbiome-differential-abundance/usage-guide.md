# Differential Abundance - Usage Guide

## Overview

Differential abundance (DA) testing identifies which individual taxa differ between groups on an amplicon ASV/feature table while respecting the compositional, relative nature of the counts. The headline workflow is a CONSENSUS of two or more compositionally-aware tools: which taxa come out "significant" depends more on the DA tool than on the biology (Nearing 2022, across 38 datasets), so a single tool's hit list is not a defensible deliverable. Run at least two of ALDEx2, ANCOM-BC2, MaAsLin2/3, LinDA, or ZicoSeq, report the intersection as high-confidence and the union as exploratory, name every tool, and disclose disagreement. A relative-abundance increase is also not an absolute increase without an external load anchor (spike-in / flow cytometry / qPCR).

This skill owns per-taxon DA on an amplicon table. Whole-community alpha/beta diversity and PERMANOVA live in diversity-analysis; the same DA math on shotgun profiler tables lives in metagenomics/metagenome-visualization; the shared compositional/closure/CLR/zero theory lives in metagenomics/abundance-estimation.

## Prerequisites

A feature table of integer counts (typically a phyloseq object) plus sample metadata. See SKILL.md's Installation section for package installs and its Tool Taxonomy for per-tool input requirements (orientation, counts vs. proportions).

## Quick Start

Tell your AI agent what you want to do:
- "Find differentially abundant taxa between treatment and control with two methods and give me the consensus"
- "Run ALDEx2 and report effect sizes with BH-adjusted q-values"
- "Run ANCOM-BC2 with age and sex as covariates and only keep hits that pass the sensitivity analysis"
- "Analyze a longitudinal study with a subject random effect using LinDA or MaAsLin2"
- "Filter taxa present in fewer than 10% of samples before testing and check the result is not sensitive to that cutoff"

## Example Prompts

### Consensus across tools
> "I have an ASV table and metadata with two treatment groups. Run at least two compositionally-aware DA methods, report the intersection as high-confidence and the union as exploratory, and tell me which tools found each taxon."

### Conservative two-group test
> "Run ALDEx2 to compare taxon abundance between healthy and diseased samples, gate on both BH-adjusted q below 0.05 and an effect-size floor, and explain the effect-size metric."

### Covariates and sensitivity
> "Run ANCOM-BC2 with age and sex as covariates, set the p-adjust method to BH, and only report hits that pass the pseudo-count sensitivity analysis (passed_ss)."

### Repeated measures
> "My samples are repeated within subjects over time. Use a DA method with a subject random effect so I do not pseudo-replicate, and cross-check the mixed-model hits across tools."

### Relative vs absolute
> "Is this taxon's increase relative or absolute? I do not have load data - explain what I can and cannot claim, and what a spike-in or qPCR anchor would add."

## Tips

The workflow steps, gotchas and thresholds above (consensus panel, effect-size + q gating, ANCOM-BC2 Holm-vs-BH, prevalence-filter sensitivity, relative-vs-absolute, DESeq2/edgeR caveat, uncorrected-test pitfall, no-settled-best-tool, organelle/contaminant filtering) are documented once, in SKILL.md - see its Decision Tree, Per-Method Failure Modes, Common Errors and Quantitative Thresholds tables.

## Related Skills

- diversity-analysis - Whole-community alpha/beta/PERMANOVA; answer "do the communities differ" first
- taxonomy-assignment - Collapse ASVs to genus/species before per-taxon testing
- amplicon-processing - Produces the ASV feature table tested here
- qiime2-workflow - The qiime composition ancombc CLI route
- metagenomics/abundance-estimation - Shared compositional/closure/CLR/zero/load-anchor theory
- metagenomics/metagenome-visualization - The same DA mechanics on shotgun profiler tables
- experimental-design/multiple-testing - FDR control and multiplicity across taxa
