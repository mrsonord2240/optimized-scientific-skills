# Permutation importance and linear models

No script. Both are short enough to show in full.

Permutation importance breaks under correlation the same way SHAP does: shuffling gene A while correlated gene B is intact lets the model recover the signal through B, so both look unimportant. Unrestricted permutation also forces predictions on points that cannot occur.

```python
from sklearn.inspection import permutation_importance
# Cluster correlated features first (see routes/shap-tree.md module map) and keep one per cluster,
# then score on HELD-OUT data, never the training rows.
r = permutation_importance(model, X_test[kept_columns], y_test, scoring='roc_auc', n_repeats=20, random_state=0)
```

- sklearn has no `conditional=` flag. For conditional permutation under correlation use R `party::cforest` with `varimp(conditional = TRUE)`.
- Always evaluate on held-out data.
- Penalized linear model: read the coefficients, the model is its own explanation. For SHAP use `shap.LinearExplainer` and choose `interventional` or `correlation_dependent` explicitly, since they give different values.
- High-stakes clinical decisions: prefer an inherently interpretable model (sparse linear, rule list), which avoids the conditional-versus-marginal ambiguity.

Done when importance is reported on held-out data over decorrelated features, with the correlation handling stated.
