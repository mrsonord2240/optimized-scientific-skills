# Perturb-seq Analysis - Usage Guide

## Overview

Perturb-seq and CROP-seq read out a pooled CRISPR screen with single-cell transcriptomes, linking a recoverable guide identity (genotype) to a transcriptome-wide phenotype. This skill treats the hard decisions: guide assignment as a mixture problem, escaper removal with Mixscape, calibrated testing with SCEPTRE, effect size with E-distance, separating compositional shifts from within-state expression, and a clear-eyed view that perturbation-prediction foundation models do not yet beat simple baselines.

## Prerequisites

See SKILL.md's Prerequisites section for exact install commands and version caveats
(`pertpy[jax]` is required, not optional, for the default guide-assignment method; `sceptre`
is installed from GitHub).

## Quick Start

Tell your AI agent what you want to do:
- "Assign guides with a mixture model, not a flat threshold"
- "Use Mixscape to remove cells that received a guide but were not perturbed"
- "Test each perturbation with SCEPTRE so false positives are calibrated"
- "Rank my perturbations by E-distance and run the E-test"
- "Tell me whether this perturbation moves cells across states or changes a state"
- "Is this foundation-model prediction actually better than a mean baseline?"

## Example Prompts

### Guide Assignment
> "Fit a Poisson-Gaussian mixture to my guide counts and assign by posterior, then show the NT contamination floor"
> "Gate doublets before treating multi-guide cells as combinatorial"

### Effective Perturbation
> "Run Mixscape to classify KO vs non-perturbed cells and report the perturbed fraction per target"
> "This target is all non-perturbed; is the gene non-functional or did the guide fail to edit?"

### Testing and Effect Size
> "Run SCEPTRE with a calibration check before the discovery analysis"
> "Compute pairwise E-distances in PCA space and run the permutation E-test against NT"
> "Do pseudobulk DESeq2 per replicate, summing raw counts, for the within-state program change"

### Composition vs Expression
> "Run Milo to test whether the perturbation shifts cell-state proportions"
> "Separate the compositional shift from the within-state expression change"

### Foundation Models
> "Benchmark this perturbation predictor on held-out whole perturbations against an additive baseline, scored on DE genes"

## What the Agent Will Do

1. Assign guides with a per-guide background/foreground mixture and report the perturbed fraction.
2. Gate doublets and decide the MOI regime (low for single-gene attribution, high for combinatorial).
3. Run Mixscape to remove non-perturbed escaper cells before any DE.
4. Choose a calibrated test (SCEPTRE conditional resampling) over naive Wilcoxon/NB.
5. Quantify effect size with E-distance in a pinned PCA embedding and the permutation E-test.
6. Run pseudobulk DE per biological replicate (summing raw counts) for the within-state program.
7. Run a differential-abundance test (Milo/scCODA) and report composition separately from expression.

Decision guidance (mixture vs threshold, which test, composition vs expression) and the
reasoning tips behind each are in SKILL.md's Governing Principle, Method Decision Tables, and
Common Errors sections — that's what the agent loads and acts on, so it's kept in one place.

## Related Skills

single-cell/preprocessing - scRNA-seq QC and normalization upstream of the screen
single-cell/doublet-detection - gating doublets before multi-guide analysis
single-cell/markers-annotation - interpreting per-perturbation DE genes
single-cell/batch-integration - multi-sample/replicate integration
crispr-screens/mageck-analysis - bulk CRISPR screen analysis (MAGeCK RRA/MLE)
crispr-screens/perturb-seq-analysis - related single-cell CRISPR screen workflow
differential-expression/deseq2-basics - pseudobulk DESeq2 testing on summed counts
pathway-analysis/go-enrichment - pathway interpretation of perturbation signatures
