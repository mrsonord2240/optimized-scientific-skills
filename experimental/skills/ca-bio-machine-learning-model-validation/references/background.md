# Background

Not needed to run an analysis. This explains choices the routes and scripts already make.

## Two root causes of inflated results

Almost every inflated result in machine learning for biology comes from one of two things: information from the test data leaked into model construction, or the same data was used both to choose and to grade a decision. A clean train/test split is necessary but not sufficient; the leak has usually happened earlier (a scaler fit on all data, a duplicate patient, ComBat across the split). The bias is largest when the true signal is weakest, which is the omics regime (p much larger than n). Selecting features on all samples can manufacture near-perfect cross-validation from pure noise.

Discrimination (AUC) and calibration are orthogonal. AUC is unchanged by any monotone transform of the score, so it cannot show whether probabilities match observed frequencies. For any decision that uses the probability itself, calibration is the property that matters.

## Choosing the scheme

| Question the estimate answers | Scheme | Why |
|-------------------------------|--------|-----|
| A new sample like the training samples, nothing tuned | Repeated stratified k-fold with the whole `Pipeline` inside | Everything is refit per fold; nothing is left for an inner loop to protect |
| Same, with tuning or model choice | Nested CV | Tuning and grading on one CV is optimistic; the optimism grows with more configurations, weaker signal, smaller n |
| A new patient | `GroupKFold` / `StratifiedGroupKFold` by patient | The unit of independence is not the row |
| A new site | Leave-one-site-out (internal-external CV) | Approximates external validation |
| Later samples | `TimeSeriesSplit`, forward-chaining | Random folds leak the future |
| Small n | `RepeatedStratifiedKFold` 5 x 10, with spread | A single CV is one high-variance draw |

Leave-one-out is high-variance and gives no AUC on a size-1 fold. The .632+ bootstrap is optimistic for learners with zero apparent error; for a single fixed model, bootstrap optimism correction is the cleaner internal validation. `StratifiedGroupKFold` is best-effort and can return single-class test folds when groups are few, which is why `scripts/cv_auc.py` stratifies on the per-patient label.

## Calibration notes

Equal-mass bins (`strategy='quantile'`) behave better than equal-width bins under imbalance. Equal-width ECE is biased: it reports error even for perfectly calibrated models. Platt scaling suits small calibration sets; isotonic needs hundreds of points. Decision curve analysis assumes calibrated probabilities; a model can have high AUC and zero net benefit at every plausible threshold.

`CalibratedClassifierCV(cv='prefit')` was deprecated in scikit-learn 1.6 and removed in 1.8; wrap the fitted model in `sklearn.frozen.FrozenEstimator`.

## Imbalance corrections

Oversampling or SMOTE changes training prevalence and shifts minority-class probabilities toward the minority; the shift grows with the resampling ratio and matters most under severe imbalance. At mild imbalance AUC and proper scores can move either way.

## Sample size and reporting

The "10 events per variable" heuristic is obsolete. The Riley minimum-sample-size framework sizes for shrinkage of at least 0.9 and precise risk estimation. TRIPOD+AI (2024) supersedes TRIPOD 2015 and asks for data-splitting and leakage controls, calibration, subgroup performance and uncertainty. PROBAST+AI is the companion risk-of-bias tool.

## Leakage by effect size (indicative)

A single global scaler may move the estimate little. Fitted quantile normalisation, ComBat/SVA, PCA or kNN/MICE imputation can move it a great deal. Feature selection on all samples is the most severe common case.
