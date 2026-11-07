# Choosing an explanation method

No script. Pick the row, then run the matching route.

| Method | Estimates | Under correlation | Use |
|--------|-----------|-------------------|-----|
| TreeSHAP `tree_path_dependent` | Conditional Shapley via tree coverage | Shifts credit among correlated split features | Fast cohort summary when conditional semantics are acceptable |
| TreeSHAP `interventional` | Marginal Shapley, features replaced from a background | Credit split can differ; cost scales with background (about 100-1000) | Default for "what does the model use", with a stated baseline |
| KernelSHAP | Model-agnostic Shapley by masking | Off-manifold masking, corrupted | Last resort for non-tree, non-net models; expensive |
| DeepSHAP, GradientSHAP | SHAP for networks against a background | Background dependent | Neural omics models |
| LinearExplainer | Exact Shapley for linear models | `interventional` and `correlation_dependent` differ | Penalized linear models, or read coefficients |
| LIME | Local sparse linear surrogate | Off-manifold, seed unstable | One local prediction only |
| Permutation importance | Score drop when a feature is shuffled | Zeroes correlated pairs | Screen on decorrelated features |
| Conditional permutation | Importance within correlated strata | Fairer | Random forest under correlation (R) |

| Scenario | Route |
|----------|-------|
| Which genes the model keys on, debugging shortcuts | `routes/shap-tree.md` with `--mode interventional` |
| Which genes are informative here (descriptive) | `routes/shap-tree.md` with `--mode tree_path_dependent`, never called reliance |
| Ranking across correlated genes | Aggregate over modules first (`routes/shap-tree.md`) |
| One local prediction | `routes/lime-local.md` or a SHAP waterfall |
| Linear model, or permutation importance | `routes/permutation-linear.md` |
| Biomarker panel | `bio-machine-learning-biomarker-discovery` |

Every method either extrapolates off the data manifold (interventional SHAP, unrestricted permutation) or leaks credit through correlation (conditional SHAP). Decide which the question tolerates.

Done when the method, its estimand and its correlation behaviour are stated with the result.
