# Regularized logistic baseline

```bash
python scripts/logistic_regression.py data.csv --label label --positive AML --out baseline_results.tsv
```

Fits scaled elastic-net logistic regression (saga, `l1_ratios` 0.1, 0.5, 0.9, 20 values of C, 5-fold CV on log loss) on a stratified 70/30 split. Writes metric/value rows: non-zero coefficients, C, l1_ratio, test AUC, Brier, calibration slope, plus `test_predictions.csv` (label, prob) for `ca-bio-machine-learning-model-validation`'s calibration route. Options: `--test-size`, `--seed`, `--l1-ratios`.

- The log-loss fit is dense, not a small signature: on Golub (2,000 probes) it kept 1,059-2,000 coefficients at the same test AUC as sparser fits. When a short list is the deliverable, add `--signature`: a lasso at the strongest penalty within one standard error of the best CV score (Golub: 45, 142 and 60 non-zero over three splits at unchanged AUC). The sparser fit shrinks probabilities and had a slightly worse Brier, so use the dense fit for probabilities and the signature for a list.
- Standardize inside the Pipeline; the penalty is scale sensitive. `scoring` is explicit, and `penalty=` is not used (removed in scikit-learn 1.10; `l1_ratio` 1 is L1, 0 is L2).
- The strongest-coefficient list is not a validated biomarker set; correlated features make selection unstable.
- Report the hold-out AUC as a fit check only. For an estimate, use `ca-bio-machine-learning-model-validation`.

Done when `baseline_results.tsv` exists and the answer states the model, non-zero count and test AUC with Brier.
