---
name: bio-machine-learning-omics-classifiers
category: Data Analysis
description: Use when building a classifier from expression, methylation, or variant data, choosing an algorithm for high-dimensional small-n data, or diagnosing a suspiciously perfect AUC.
tool_type: python
primary_tool: sklearn
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: numpy 2.5, pandas 3.0, scikit-learn 1.9, xgboost 3.4, imbalanced-learn 0.14 (the bundled scripts were run against these).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

Three high-risk drifts: XGBoost moved `early_stopping_rounds` from `fit()` to the constructor (deprecated 1.6, removed from `fit()` in 2.1); scikit-learn 1.8 deprecated `penalty=` on `LogisticRegression`/`LogisticRegressionCV` (removed in 1.10; leave it unset and use `l1_ratio`/`l1_ratios` + `C`: `l1_ratio=1` is L1, `0` is L2, in between is elastic net) and made `LogisticRegressionCV` warn until `scoring=` and `use_legacy_attributes=` are set explicitly; `CalibratedClassifierCV(cv='prefit')` was deprecated in 1.6 and removed in 1.8 (use `FrozenEstimator`). If code throws TypeError/FutureWarning, switch to the constructor / `l1_ratio` / `FrozenEstimator` form.

# Classification Models for Omics Data

**"Build a classifier from my expression data"** -> Start with a regularized linear model (often the ceiling in p>>n), check for batch shortcuts, and treat the probability -- not the label -- as the product.
- Linear (often best): `LogisticRegression(solver='saga', l1_ratio=0.5)` (elastic net)
- Trees when nonlinear/interaction signal: `RandomForestClassifier`, `xgboost.XGBClassifier`
- Imbalance: keep the natural class prior and tune the threshold; class weights and SMOTE shift predicted risk, so use them only for hard-label problems and recalibrate if a probability is needed

## The Single Most Important Modern Insight -- In p>>n, Simple Often Wins and the Probability, Not the Label, Is the Product

Omics classification almost always lives in p>>n (thousands of features, tens-to-hundreds of samples). Two counterintuitive consequences follow. First, more flexible is not better: with n in the dozens the variance of a flexible learner dominates, the full covariance is singular so QDA/full-LDA are undefined, and simple diagonal/linear methods match or beat elaborate ones (Dudoit 2002). "Random forest is the obvious choice for expression data" is a myth -- SVM/regularized logistic frequently win on microarray-style problems (Statnikov 2008), and gradient-boosted trees beat deep nets on tabular/omics data (Grinsztajn 2022). Regularization is the load-bearing wall, not a tuning nicety.

Second, in diagnostic/prognostic use the probability is the product, not the label -- which makes calibration, not accuracy, the thing that breaks silently (Van Calster 2019). A model can rank perfectly (AUC 0.9) and still output dishonest risks. And the most common cause of a beautiful AUC is not skill but a batch artifact: if batch correlates with the outcome, the classifier learns the cleaner technical signal and the performance collapses on any independent cohort.

## Algorithm Choice for p>>n

Not executed here: the LightGBM, CatBoost, linear SVM and DLDA/nearest-centroid guidance (this table, the Platt-scaling advice, the preprocessing notes) is literature-based; no script or run covers it. Everything about elastic-net logistic, RF, XGBoost, calibration and batch checks was executed.

| Model | Wins when | Overfits / fails when | Calibration | Scaling |
|-------|-----------|------------------------|-------------|---------|
| L1 logistic (lasso) | Sparse signal, want a small signature | Correlated features -> unstable selection; >n true signals | Good (proper loss); shrinks toward base rate | Standardize |
| L2 / elastic-net logistic | Many small correlated effects; omics default | Needs C (and l1_ratio) tuning | Good; preferred when calibration matters | Standardize |
| DLDA / nearest-centroid | Tiny n, roughly linear (Dudoit 2002) | Strong interactions; non-Gaussian | Crude; recalibrate | Variance-scaled |
| Linear SVM | High-dim linear separability (Statnikov 2008) | Heavy overlap; needs C | No native probabilities -- Platt-scale `decision_function` | Critical |
| Random forest | Nonlinear/interaction signal; robust baseline | Sparse-linear signal; tiny n; OOB-as-test leakage | Bagged votes bounded away from 0 and 1 | Scale-invariant |
| GBDT (XGBoost/LightGBM) | Best general tabular performer | Tiny n + deep/many rounds; needs early stopping | Log-loss overfitting tends to overconfident extremes | Scale-invariant |
| Tabular deep nets | Very large n; multimodal/transfer | Typical omics n -> loses to GBDT | Variable; often needs temperature scaling | Standardize |

Tree ensembles are often miscalibrated and the direction depends on the learner and loss: classic boosted ensembles push probabilities toward 0.5 (sigmoid distortion; Niculescu-Mizil 2005), bagged forests are comparatively well-calibrated but bound their votes away from 0 and 1, and modern gradient boosting trained to log-loss for many rounds tends to overfit toward overconfident extremes. Check a reliability curve and recalibrate rather than assuming a direction.

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Default omics classifier, want a signature | Elastic-net logistic | Often the ceiling in p>>n; sparse + grouping; well-calibrated |
| Suspected nonlinear/interaction (epistasis, thresholds) | Random forest then XGBoost with early stopping | Trees capture interactions; benchmark vs the linear baseline |
| Probabilities will drive a clinical decision | Linear model + calibration check; recalibrate if needed | Probability is the product; AUC is blind to calibration |
| Class imbalance, probability is the output | No reweighting or resampling; tune the threshold | Class weights and SMOTE shift predicted risk (measured below); ranking changes only for tree ensembles |
| Class imbalance, hard label only | Class weights or SMOTE inside the training fold | Sensitivity matters more than the probability scale; recalibrate if a probability is later needed |
| Mixed continuous + categorical features | `ColumnTransformer` (scale continuous, encode categorical) | Different feature types need different handling |
| Missing values, especially below-detection | XGBoost/LightGBM native NaN handling | Missingness is often informative (MNAR) |
| Considering a deep net | Only at very large n or multimodal/raw inputs | GBDT beats deep on engineered omics matrices (Grinsztajn 2022) |
| Need unbiased performance / nested CV / calibration metrics | -> machine-learning/model-validation | Evaluation is its own discipline |
| Time-to-event outcome | -> machine-learning/survival-analysis | Censoring needs survival models, not classifiers |

## Core Workflow: Regularized Logistic First

**Goal:** A calibrated, interpretable baseline that is often the best omics classifier.

**Approach:** Standardize inside a Pipeline and fit elastic-net logistic with a cross-validated penalty; the L2 component keeps correlated genes together. Name `scoring` explicitly: it decides both the probability quality and how many coefficients survive.

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV

# saga supports mixed L1/L2 (0 < l1_ratio < 1); C = 1/lambda (small C = strong shrinkage). Standardize: the penalty is scale-sensitive.
# neg_log_loss selects the penalty that gives the best probabilities, which is a weak penalty (see below).
clf = LogisticRegressionCV(solver='saga', l1_ratios=[0.1, 0.5, 0.9], Cs=20, cv=5, max_iter=10000,
                           scoring='neg_log_loss', use_legacy_attributes=False, random_state=0)
pipe = Pipeline([('scaler', StandardScaler()), ('clf', clf)])
pipe.fit(X_train, y_train)
```

**This fit is not a small signature.** Measured (70/30 Golub splits, 2,000 highest-variance probes; synthetic 500 features, 15 true): `neg_log_loss` keeps 1,059-2,000 of 2,000 Golub coefficients and 70-82 of 500 synthetic ones (13-15 of the 15 true found); `accuracy` keeps 215-481 and `roc_auc` 157-228 of 2,000 on Golub. Test AUC is the same either way. When a short list is the deliverable, fit a lasso and take the strongest penalty within one standard error of the best CV score:

```python
import numpy as np
from sklearn.linear_model import LogisticRegression

lasso = LogisticRegressionCV(solver='saga', l1_ratios=[1.0], Cs=20, cv=5, max_iter=10000,
                             scoring='neg_log_loss', use_legacy_attributes=False, random_state=0)
Pipeline([('scaler', StandardScaler()), ('clf', lasso)]).fit(X_train, y_train)
folds = lasso.scores_[:, 0, :]                                  # (cv folds, Cs)
mean, se = folds.mean(0), folds.std(0, ddof=1) / np.sqrt(folds.shape[0])
c_1se = lasso.Cs_[np.flatnonzero(mean >= mean.max() - se[mean.argmax()]).min()]   # Cs_ ascends: smallest C = strongest penalty
signature = Pipeline([('scaler', StandardScaler()),
                      ('clf', LogisticRegression(solver='saga', l1_ratio=1.0, C=c_1se, max_iter=10000, random_state=0))]).fit(X_train, y_train)
print(int((signature['clf'].coef_ != 0).sum()), 'non-zero coefficients')
```

Measured on the same data: Golub 45, 142 and 60 non-zero coefficients over three splits at unchanged test AUC (1.00, 0.99, 1.00); synthetic 26-85 non-zero with 13-14 of 15 true found, AUC within 0.012, Brier 0.006-0.009 worse in all three seeds (the stronger penalty shrinks probabilities). Choose the dense fit for probabilities and this one for a list; Golub is near-separable, so it shows sparsity only, not discrimination.

## Tree Ensembles

**Goal:** Capture nonlinear and interaction structure when the linear baseline leaves signal on the table.

**Approach:** Random forest needs no scaling and is a robust baseline; XGBoost needs a low learning rate, shallow depth, and a rule for the number of rounds (early stopping, set in the constructor in 2.x) to avoid overfitting tiny n. Keep three disjoint parts: train, an early-stopping validation set, and a test set that only the final score touches. Early-stopping data is model-selection data: never report its metric as performance.

```python
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier

X_dev, X_te, y_dev, y_te = train_test_split(X, y, test_size=0.3, stratify=y, random_state=0)
X_tr, X_val, y_tr, y_val = train_test_split(X_dev, y_dev, test_size=0.3, stratify=y_dev, random_state=0)

rf = RandomForestClassifier(n_estimators=500, max_features='sqrt', min_samples_leaf=3,
                            n_jobs=-1, random_state=0)

# XGBoost 2.x: early_stopping_rounds and eval_metric go in the CONSTRUCTOR, not fit().
# scale_pos_weight is omitted on purpose: like class_weight it reweights the prior and
# shifts predicted risk -- use it only for hard-label problems (see Class Imbalance).
xgb = XGBClassifier(n_estimators=2000, learning_rate=0.03, max_depth=4, subsample=0.8,
                    colsample_bytree=0.5, reg_lambda=1.0,
                    early_stopping_rounds=50, eval_metric='aucpr', n_jobs=-1, random_state=0)
xgb.fit(X_tr, y_tr, eval_set=[(X_val, y_val)], verbose=False)      # NaN handled natively (missing=np.nan)
print('stopped at round', xgb.best_iteration)                     # report scores on (X_te, y_te) only
```

A small validation set stops early and flatters itself: with 120 training samples and a 15-sample validation set the median best round was 4.5 and validation AUCPR averaged 0.83 against a test AUC of 0.61 (60 samples: 0.72 vs 0.65; Golub with 15: validation AUCPR 1.0 in every split, best round 4-18). At small n, choose the rounds by cross-validation on the development data instead (measured at 120 training samples: test AUC 0.69 vs 0.61-0.66 for early stopping on 15-120 validation samples):

```python
from sklearn.model_selection import GridSearchCV

search = GridSearchCV(XGBClassifier(learning_rate=0.03, max_depth=4, subsample=0.8, colsample_bytree=0.5,
                                    n_jobs=-1, random_state=0),
                      {'n_estimators': [100, 200, 400]}, cv=5, scoring='roc_auc').fit(X_dev, y_dev)
```

`python scripts/rf_xgboost_classifier.py` compares elastic-net logistic, RF, and XGBoost on synthetic p>>n data (600 samples, 1,500 features; XGBoost stops at round 223 and reaches test AUC 0.83) by test-split AUC and Brier and prints each model's probability spread. At 300 samples the same XGBoost setting stops at round 0 (AUC 0.60, a null model); the script prints a warning then.

## Detecting Batch Shortcut Learning

**Goal:** Rule out that a high AUC is a batch artifact rather than biology.

**Approach:** Try to predict the batch from the features, then hold out whole batches; if batch is confounded with the outcome, no correction rescues the design (Soneson 2014) -- fix it at the design stage. `scripts/batch_checks.py` works for any number of batches from 2.

```python
import sys
sys.path.insert(0, 'scripts')       # relative to the cwd: run from the Skill directory, or use the absolute path of its scripts/
import pandas as pd
from scipy.stats import chi2_contingency
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from batch_checks import batch_predictability, leave_one_batch_out

probe = Pipeline([('scaler', StandardScaler()), ('clf', LogisticRegression(max_iter=5000))])

# 1. Can the features predict the BATCH (one-vs-rest AUC)? If yes, batch is a strong axis and the label model is suspect.
print(f'Batch predictability AUC: {batch_predictability(probe, X, batch_labels):.2f} (high = shortcut risk)')

# 2. Is the outcome associated with batch by design?
print('label vs batch p:', chi2_contingency(pd.crosstab(y, batch_labels))[1])

# 3. Leave-one-batch-out: each batch is predicted by a model that never saw it.
pooled, per_batch = leave_one_batch_out(probe, X, y, batch_labels)
print(per_batch.to_string(index=False))
print(f'Leave-one-batch-out mean per-batch AUC: {per_batch.auc.mean():.2f} (pooled {pooled:.2f})')
```

Read the per-batch mean, not the pooled AUC: pooling also ranks samples across batches, so a batch-level score shift that tracks each batch's case rate still scores. A held-out batch with one class has no AUC and shows as NaN with a note (its predictions still enter the pooled number); a fold whose training batches hold one class is skipped and listed. Measured on synthetic zero-signal data confounded with batch (random-split label AUC 0.82 / 0.88 / 0.86 for 2 / 3 / 6 batches): leave-one-batch-out mean per-batch AUC 0.54 / 0.70 / 0.45, so the estimate is noisy at any batch count (the spread of a single leave-one-batch-out mean was 0.09-0.10 at 2, 3 and 6 batches), and one run is not a verdict; batch predictability was 1.00 in all three.

`python scripts/logistic_regression.py` fits an elastic-net logistic model and demonstrates the batch-confounding artifact on synthetic data.

## Class Imbalance: What Works and What Fails

**Goal:** Handle a rare positive class without destroying the probabilities.

**Approach:** When the probability is the output, change nothing about the training prior: no `class_weight`, no `scale_pos_weight`, no SMOTE. Pick the operating threshold by cost on a validation fold; the sensitivity that resampling buys is recoverable that way (van den Goorbergh 2022). Reweighting and resampling always change the probability scale; for a logistic model that is all they change, while for a tree ensemble they can change the ranking too (better or worse, so test it). Either way recalibrate before reading a risk. If a hard-label problem justifies them, use an `imblearn` Pipeline so only training folds are resampled, and recalibrate on held-out data at the natural prevalence before reading any output as a risk.

Measured on a known logistic model at 8% prevalence (`python scripts/calibration_check.py`; 1,500 train, 1,000 calibration, 20,000 test, mean of 5 seeds): unweighted logistic predicts 0.99x the observed prevalence (calibration slope 0.93, AUC 0.77); `class_weight='balanced'` predicts 4.2x (slope 0.72, AUC 0.77, Brier 0.177 vs 0.070); the balanced model recalibrated on the 1,000 held-out samples predicts 1.11x (Brier 0.072). A random forest (300 features, 500 train) predicted 2.6x the prevalence with `class_weight='balanced'` and 1.16x without. The 1,000 calibration samples held about 80 events. Ranking by family (own generator, 8% prevalence, 6 seeds, test AUC plain vs balanced): logistic 0.826 vs 0.820, so the scale moves and the ranking does not; random forest 0.776 vs 0.812, higher in 6 of 6 seeds, so the weights changed the ranking (in this setup for the better). Balanced weights still left RF risks at 2.4x the prevalence.

```python
from imblearn.pipeline import Pipeline as ImbPipeline   # NOT sklearn's Pipeline
from imblearn.over_sampling import SMOTE
from sklearn.linear_model import LogisticRegression

# Hard-label problems only. SMOTE's fit_resample runs during fit on the train fold, no-op on transform.
imb = ImbPipeline([('smote', SMOTE(random_state=0)),
                   ('clf', LogisticRegression(max_iter=5000))])
# Risk models: fit the plain classifier and choose the threshold by cost on a validation fold.
```

## Probability Calibration

**Goal:** Ensure a "0.9" means a 90% risk, not just a high rank.

**Approach:** Tree ensembles are often miscalibrated in a learner-dependent direction (see the note under the algorithm table), so check a reliability curve and recalibrate on a disjoint fold. Logistic regression optimizes a proper scoring rule and is usually best-calibrated out of the box. See machine-learning/model-validation for reliability curves, Brier, and the full protocol.

Method by calibration-set size: use `sigmoid` (Platt) by default; use `isotonic` only with at least 1,000 calibration samples and 100 in the rarer class. Measured on RF and XGBoost bases (50 features, 10% and 50% prevalence, 40 repeats, 10,000 test samples): isotonic had the lower test Brier in at most 25% of repeats for RF at any size up to 1,000, and in at most 50% for XGBoost up to 200 samples (RF at 50% prevalence: Brier 0.2545 isotonic vs 0.2326 sigmoid at 20 samples, 0.2308 vs 0.2221 at 100). At 1,000 samples isotonic was within 0.002 of sigmoid for RF and won 55-80% of repeats for XGBoost (93-100% at 3,000). When the base model separates the calibration set perfectly, isotonic returns only 0 and 1 (Golub, 20 samples: two distinct values, Brier 0.000; reproduced in `scripts/calibration_check.py`). That Brier describes the calibrator's own fit, not performance: never report it, and treat three or fewer distinct output values as a failed calibration. A small calibration set can also lose to the raw probabilities (RF, 30% prevalence, 20 samples: Brier 0.1985 raw vs 0.2089 sigmoid vs 0.2356 isotonic; at 100 samples sigmoid only ties raw at 0.1980), so compare against the uncalibrated Brier and recalibrate only when it improves.

```python
import sys
sys.path.insert(0, 'scripts')                            # same path rule as above
from calibration_check import recalibrate, calibration_table

# recalibrate wraps FrozenEstimator (sklearn >=1.6; cv='prefit' deprecated), picks sigmoid/isotonic by the rule
# above, and warns if the fitted calibrator has <= 3 distinct probabilities.
calibrated = recalibrate(rf.fit(X_tr, y_tr), X_cal, y_cal)    # X_cal disjoint from train and test
p_te = calibrated.predict_proba(X_te)[:, 1]
print(calibration_table(y_te, p_te))                          # mean predicted vs observed risk per bin; report this and the Brier on X_te only
```

## Preprocessing, Scaling, Encoding, Missing Data

- Inside the CV fold (refit per fold via Pipeline): feature selection, scaling, log/VST/quantile normalization, ComBat batch correction, imputation, PCA. Fitting any of these on the full matrix is leakage (machine-learning/model-validation).
- Scale-sensitive: regularized logistic, SVM, kNN, PCA, neural nets -- standardize. Scale-invariant: trees/RF/GBDT.
- Variant features: additive 0/1/2 ordinal for trees/additive logistic; one-hot for non-additive (dominant/recessive) effects; CatBoost-style encoding for high-cardinality (HLA).
- Missing data: XGBoost/LightGBM learn a default split direction for NaN; in omics missingness is often MNAR (below detection limit) and a missingness indicator can carry signal -- naive zero-fill conflates "absent" with "not measured."

## Hyperparameters That Matter

| Model | Tune these | Leave default |
|-------|-----------|---------------|
| Logistic | `C` (log-spaced), `l1_ratio` (`class_weight` only for hard-label problems) | solver (saga for elasticnet) |
| Random forest | `max_features`, `min_samples_leaf`, `max_depth` (cap for tiny n) | `n_estimators` (more is safe; 500-1000) |
| XGBoost | `learning_rate`+`n_estimators`+early stopping, `max_depth` (3-6), `subsample`, `colsample_bytree`, `reg_lambda` | most others |

## Failure Modes

Read [`references/failure-modes.md`](references/failure-modes.md) when a classifier AUC looks too good, probabilities look dishonest, or imbalance handling is in question.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Try a regularized linear model first | Dudoit 2002; Statnikov 2008 | Simple often beats complex in p>>n |
| GBDT over deep nets for tabular omics | Grinsztajn 2022 | Trees handle uninformative features and non-rotational data |
| Do not reweight or resample for risk models | van den Goorbergh 2022; Carriero 2025 | Both shift predicted risk; logistic AUC unchanged, RF AUC changes |
| XGBoost: low LR + many rounds + early stopping | field standard | Prevents overfitting tiny n |
| Report AUPRC + MCC under imbalance | Saito 2015; Chicco 2020 | Accuracy and ROC-AUC mislead when positives are rare |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| XGBoost `early_stopping_rounds` TypeError in `fit()` | Moved to constructor in 2.x | Pass it (and `eval_metric`) in `XGBClassifier(...)` |
| `penalty=` FutureWarning | Deprecated in sklearn 1.8, removed in 1.10 | Leave `penalty` unset; use `l1_ratio` (`l1_ratios` in the CV class) + `C` |
| Mixed L1/L2 solver error | Only saga supports `0 < l1_ratio < 1` | `solver='saga'` + `l1_ratio` |
| `CalibratedClassifierCV(cv='prefit')` raises | Deprecated in 1.6, removed in 1.8 | Wrap in `FrozenEstimator` |
| 95% accuracy but useless model | Imbalance + accuracy metric | Report AUPRC/MCC; check the confusion matrix |

## References

- Dudoit S, Fridlyand J, Speed TP. 2002. Comparison of discrimination methods for the classification of tumors using gene expression data. *J Am Stat Assoc* 97:77-87.
- Chawla NV, Bowyer KW, Hall LO, Kegelmeyer WP. 2002. SMOTE: Synthetic Minority Over-sampling Technique. *J Artif Intell Res* 16:321-357.
- Zou H, Hastie T. 2005. Regularization and variable selection via the elastic net. *J R Stat Soc B* 67:301-320.
- Niculescu-Mizil A, Caruana R. 2005. Predicting good probabilities with supervised learning. *Proc 22nd ICML* 625-632.
- Diaz-Uriarte R, Alvarez de Andres S. 2006. Gene selection and classification of microarray data using random forest. *BMC Bioinformatics* 7:3.
- Statnikov A, Wang L, Aliferis CF. 2008. A comprehensive comparison of random forests and support vector machines for microarray-based cancer classification. *BMC Bioinformatics* 9:319.
- Soneson C, Gerster S, Delorenzi M. 2014. Batch effect confounding leads to strong bias in performance estimates obtained by cross-validation. *PLoS ONE* 9:e100335.
- Saito T, Rehmsmeier M. 2015. The precision-recall plot is more informative than the ROC plot when evaluating binary classifiers on imbalanced datasets. *PLoS ONE* 10:e0118432.
- Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW. 2019. Calibration: the Achilles heel of predictive analytics. *BMC Med* 17:230.
- Chicco D, Jurman G. 2020. The advantages of the Matthews correlation coefficient (MCC) over F1 score and accuracy in binary classification evaluation. *BMC Genomics* 21:6.
- Shwartz-Ziv R, Armon A. 2022. Tabular data: deep learning is not all you need. *Inf Fusion* 81:84-90.
- Grinsztajn L, Oyallon E, Varoquaux G. 2022. Why do tree-based models still outperform deep learning on typical tabular data? *NeurIPS Datasets and Benchmarks*.
- van den Goorbergh R, van Smeden M, Timmerman D, Van Calster B. 2022. The harm of class imbalance corrections for risk prediction models. *J Am Med Inform Assoc* 29:1525-1534.
- Carriero A, Luijken K, de Hond A, Moons KGM, van Calster B, van Smeden M. 2025. The harms of class imbalance corrections for machine learning based prediction models. *Stat Med* 44:e10320.

## Related Skills

- machine-learning/model-validation - Nested CV, calibration, and net benefit for the trained classifier
- machine-learning/biomarker-discovery - Select features before modeling (inside the CV fold)
- machine-learning/prediction-explanation - Interpret the classifier and detect shortcuts with SHAP
- machine-learning/survival-analysis - Time-to-event outcomes that classifiers cannot handle
- differential-expression/batch-correction - Batch correction done design-aware, not across the split
- expression-matrix/normalization - Per-sample normalization that is safe outside the CV fold
