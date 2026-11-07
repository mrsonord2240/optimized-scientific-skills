# Random forest and XGBoost against the linear baseline

```bash
python scripts/rf_xgboost_classifier.py data.csv --label label --positive AML --out trees_results.tsv
```

Compares elastic-net logistic, random forest (500 trees) and XGBoost (learning rate 0.03, depth 4, early stopping) with three disjoint splits: train, a validation split used only for early stopping, and a test split that every reported number comes from. Output: test AUC, Brier and probability spread per model.

- `best_round=0` in the note marks XGBoost a null model; do not read that row. With a validation split under about 50 samples early stopping flatters itself (validation AUCPR 0.83 against test AUC 0.61 at 120 training samples); rerun with `--xgb-rounds cv` to pick the rounds by 5-fold CV.
- XGBoost 2.x and later take `early_stopping_rounds` and `eval_metric` in the constructor, not `fit()`.
- NaN in the features: the script skips logistic and forest and runs XGBoost, which handles NaN natively. In omics, missingness is often informative (below detection); zero-fill conflates absent with not measured.
- Escalate only if trees clearly beat the linear model on the test split. A wide probability spread is not by itself miscalibration; use `routes/calibration.md`.
- Never report the validation metric as performance.

Done when `trees_results.tsv` exists and the answer says whether the trees beat the linear baseline on the test split.
