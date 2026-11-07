# Probabilities read as risks

```bash
python scripts/calibration_check.py data.csv --label label --positive X --model rf --out calibration_results.tsv
```

Same script as `routes/imbalance.md`. It prints the reliability table of the unweighted model on the test split and the recalibrated variant. For a reliability curve, Brier and slope on an existing predictions file, use `ca-bio-machine-learning-model-validation`'s calibration route.

- Logistic regression optimizes a proper scoring rule and is usually best calibrated out of the box. Tree ensembles are miscalibrated in a learner-dependent direction (forests bound votes away from 0 and 1; classic boosting is sigmoid-distorted; log-loss boosting over many rounds goes overconfident), so check, do not assume.
- Recalibrate only on a split disjoint from training and test. The script uses `CalibratedClassifierCV` around `FrozenEstimator` (`cv='prefit'` was removed in scikit-learn 1.8).
- Method: sigmoid (Platt) by default; isotonic only with at least 1,000 calibration samples and 100 in the rarer class. Isotonic beat sigmoid on Brier in at most 25% of repeats for forests up to 1,000 samples.
- A calibrator with three or fewer distinct outputs is degenerate (isotonic on 20 separable samples gives only 0 and 1, Brier 0.000). Never report its scores; the script warns.
- Compare against the uncalibrated Brier: at about 20 calibration samples recalibration can score worse than the raw probabilities.

Done when mean predicted against observed risk by bin, and the Brier on the test split, are reported from data the model and the calibrator never saw.
