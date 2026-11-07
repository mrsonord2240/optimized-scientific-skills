# Rare positives and metric choice

No script. Compute the headline metric with `routes/fixed-pipeline.md` (ROC AUC), then add the rows below that the request needs.

| Metric | Use for | Rule |
|--------|---------|------|
| Accuracy | Never as the headline under imbalance | 5% prevalence: "always negative" scores 95% |
| ROC AUC | Discrimination | Blind to calibration |
| Average precision (AUPRC) | Rare positives | State the prevalence as the baseline, not 0.5 |
| Brier, log-loss | Probabilities are used | Not comparable across prevalences without scaling |
| MCC | One balanced threshold summary | Still depends on the threshold |
| F1 | Retrieval-style problems | Ignores true negatives; assumes a cost ratio |

- Choose the operating threshold on a separate fold (or by net benefit, `routes/decision-curve.md`), then report the metric at that locked threshold once.
- Never report the best F1 or accuracy over thresholds, and never choose the threshold on the test fold.
- Prefer threshold-free curves (ROC, PR, calibration) plus one pre-specified operating point.
- SMOTE or oversampling for a risk model shifts predicted probabilities away from the true prevalence. Keep class-preserving training and pick the threshold from utility. If resampling is justified, put it inside an `imblearn.pipeline.Pipeline` (train folds only), then check calibration and recalibrate at the original prevalence.
- Use `RepeatedStratifiedKFold` so every fold holds positives.

Done when the metric, its baseline and the threshold rule are stated.
