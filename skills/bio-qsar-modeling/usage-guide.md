# QSAR Modeling Usage Guide

## Overview

Build target-specific QSAR / QSPR models from in-house assay data by comparing chemprop D-MPNN, RandomForest + ECFP4, and, when justified, transformer-based MolFormer / Uni-Mol / ChemBERTa under the same validation design. Apply OECD 5 principles with deployment-relevant splits, ensemble uncertainty, applicability-domain assessment, and conformal prediction where its assumptions and interface are satisfied.

## Quick Start

Tell the AI agent what to do:
- "Train chemprop D-MPNN classifier on this hERG dataset with scaffold split"
- "Build RandomForest QSAR on 150 compounds (small data); use ECFP4 + nested CV"
- "Add conformal prediction intervals to my chemprop classifier"
- "Apply applicability domain filter (kNN distance) to new predictions"
- "Compute SHAP atomic attribution to interpret hERG predictions"

## Example Prompts

### Small-data baseline
> "Build RF + ECFP4 binary classifier on 150-compound dataset (hERG_blocker column). 5-fold scaffold-split cross-validation. Report AUC, sensitivity, specificity, applicability domain by kNN Tanimoto."

### chemprop classifier
> "Train chemprop 2.x D-MPNN ensembles on hERG_data.csv. Use a scaffold-balanced split, five replicates, and five ensemble members per replicate. Include current `rdkit_2d` features. Report test AUC and calibration, keeping replicate and ensemble aggregation explicit."

### Conformal prediction
> "Use MAPIE with my sklearn QSAR baseline to target 90% marginal coverage, or implement and test an sklearn-compatible wrapper before applying it to chemprop. State the exchangeability assumptions."

### Multi-task CYP inhibition
> "Train chemprop multitask classifier on CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4. Scaffold split. Compare to per-target single-task models."

## What the Agent Will Do

1. Standardize input data (recommend `chemoinformatics/molecular-standardization`).
2. Choose and document a split that represents deployment (for example scaffold, time, external-series, or prospective).
3. Featurize with ECFP4 + RDKit 2D descriptors (for chemprop hybrid).
4. Train replicated ensembles for chemprop or bootstrap models for RF, keeping validation partitions independent of model selection.
5. Compute test metrics + calibration plots.
6. Build applicability domain assessment (kNN Tanimoto distance or conformal prediction).
7. Optionally compute SHAP feature importance.

## Related Skills

- chemoinformatics/molecular-descriptors - Featurization choices
- chemoinformatics/molecular-standardization - Mandatory upstream
- chemoinformatics/scaffold-analysis - Scaffold split implementation
- chemoinformatics/admet-prediction - ADMET-specific QSAR
- chemoinformatics/generative-design - QSAR as scoring component
- machine-learning/model-validation - General ML validation
- machine-learning/biomarker-discovery - Adjacent ML approaches
