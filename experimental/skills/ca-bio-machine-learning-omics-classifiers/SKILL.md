---
name: ca-bio-machine-learning-omics-classifiers
category: Data Analysis
description: Use when building a classifier from expression, methylation or variant data, choosing an algorithm for high-dimensional small-n data, or diagnosing a suspiciously perfect AUC.
tool_type: python
primary_tool: sklearn
license: MIT
author: GPTomics
---

# Classifiers for Omics Data

Pick the one row that matches your task. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your task | Read |
|-----------|------|
| Build the default classifier, or a short feature signature | `routes/baseline.md` |
| Check whether a linear model leaves signal for random forest or XGBoost | `routes/trees.md` |
| The AUC looks too good, or cases and controls were processed in different batches | `routes/batch-check.md` |
| Rare positive class | `routes/imbalance.md` |
| Predicted probabilities will be read as risks | `routes/calibration.md` |
| Which algorithm for this data, preprocessing, missing values, hyperparameters | `routes/algorithm-choice.md` |
| An unbiased performance estimate, nested CV, leakage audit, decision curves | `ca-bio-machine-learning-model-validation` |

Rules on every route:

- Input is a CSV with one row per sample, a binary label column and numeric feature columns (loaded by `scripts/omics_data.py`, which refuses a non-binary label, non-numeric columns and missing values). Name the positive class with `--positive`; name id or batch columns with `--drop`.
- In p much larger than n, start with a regularized linear model. Random forest is not the obvious choice for expression data.
- Scaling, selection, imputation and batch correction are fit on training folds only. The scripts use a Pipeline; never fit them on the full matrix first.
- When the probability is the product, do not use class weights, `scale_pos_weight` or SMOTE; tune the threshold. Reweighting shifts predicted risk.
- A single hold-out AUC from these scripts is a fit check, not the performance estimate.
- Time-to-event outcomes need `bio-machine-learning-survival-analysis`, not a classifier; feature selection is `bio-machine-learning-biomarker-discovery`; SHAP is `ca-bio-machine-learning-prediction-explanation`.

Tested with scikit-learn 1.9.1, xgboost 3.4.2, numpy 2.5.3, pandas 3.0.6, scipy 1.18.1; imbalanced-learn not installed in the checked environment.
