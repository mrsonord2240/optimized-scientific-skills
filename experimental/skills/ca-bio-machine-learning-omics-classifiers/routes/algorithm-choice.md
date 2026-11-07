# Algorithm choice and preprocessing

No script. Pick by the data, then go to the route that runs it.

| Model | Wins when | Fails when | Scaling |
|-------|-----------|------------|---------|
| L1 logistic | Sparse signal, small signature | Correlated features make selection unstable | Standardize |
| L2 or elastic-net logistic | Many small correlated effects; the omics default; calibrated | Needs C and l1_ratio tuned | Standardize |
| DLDA, nearest centroid | Tiny n, roughly linear | Strong interactions | Variance-scaled |
| Linear SVM | High-dimensional linear separability | Heavy overlap; no native probabilities (Platt-scale `decision_function`) | Critical |
| Random forest | Nonlinear or interaction signal; robust baseline | Sparse linear signal, tiny n, out-of-bag error read as a test estimate | Invariant |
| XGBoost or LightGBM | Best general tabular performer | Tiny n with deep trees or many rounds | Invariant |
| Deep nets on tables | Very large n, multimodal | Typical omics n; gradient boosting beats them | Standardize |

The SVM, DLDA, LightGBM and CatBoost rows are literature guidance and have no script.

| Situation | Route |
|-----------|-------|
| Default, or a signature wanted | `routes/baseline.md` |
| Suspected interactions or thresholds | `routes/trees.md` |
| Probability drives a decision | `routes/calibration.md` |

- Inside the CV fold, refit per fold: feature selection, scaling, log or VST or quantile normalization, ComBat, imputation, PCA. Fitting any of these on the full matrix is leakage.
- Variant features: additive 0/1/2 for trees and additive logistic; one-hot for dominant or recessive effects.
- Missing values: XGBoost and LightGBM learn a default direction for NaN.
- Tune for logistic: `C` (log-spaced), `l1_ratio`. Forest: `max_features`, `min_samples_leaf`, `max_depth`. XGBoost: learning rate with rounds and early stopping, `max_depth` 3-6, `subsample`, `colsample_bytree`, `reg_lambda`.
- Out-of-bag error after selecting features on the full data is optimistic.

Done when the algorithm is chosen with its reason and the matching route has been run.
