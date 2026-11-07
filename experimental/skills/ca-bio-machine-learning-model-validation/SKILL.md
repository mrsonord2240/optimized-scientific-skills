---
name: ca-bio-machine-learning-model-validation
category: Data Analysis
description: Use when estimating how well a classifier trained on omics or biomedical samples will perform on new samples or new patients.
tool_type: python
primary_tool: sklearn
license: MIT
author: GPTomics
---

# Honest Model Validation

Pick the one row that matches your data. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your situation | Read |
|----------------|------|
| One sample per patient, pipeline and feature count fixed in advance | `routes/fixed-pipeline.md` |
| Several samples per patient or donor | `routes/grouped.md` |
| Feature count, regularisation or model type is being chosen or tuned | `routes/nested-cv.md` |
| The model must work at a new site or on later samples | `routes/site-or-time.md` |
| Predicted probabilities will be used, not just ranks | `routes/calibration.md` |
| Does the model help a decision at some risk threshold | `routes/decision-curve.md` |
| Rare positives, or which metric to headline | `routes/imbalance-metrics.md` |
| Check an existing analysis for leakage | `routes/leakage-audit.md` |
| What to state in the write-up, sample size, external validation | `routes/write-up.md` |

Rules on every route:

- Input is a CSV with one row per sample, a binary label column, and numeric feature columns.
- Selection, scaling, imputation and batch correction are fit on training folds only. The scripts do this. Never fit them on the full matrix first.
- Split by the unit that must be independent (patient, donor, site, time), never by row when rows share one.
- Never choose and grade on the same data: no tuning, threshold or recalibration on the evaluation fold.
- Report the mean over repeats as the estimate and its SD as the uncertainty. Choose folds and repeats once, before looking. Do not add permutation tests, learning curves or more repeats unless asked.
- An AUC near 0.5 on noise-like features is the answer. Do not tune it away.

Tested with scikit-learn 1.9.1, pandas 2.3.3, numpy 2.5.3.
