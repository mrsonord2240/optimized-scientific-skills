# Leakage audit of an existing analysis

No script. Read the user's code and check each row. Fix every hit, then re-estimate with the route that fits.

| Leak | Look for | Fix |
|------|----------|-----|
| Preprocessing | `fit_transform` on all data before splitting: scaler, quantile or library normalisation, ComBat/SVA, PCA, imputation, VST | Put every transform in the `Pipeline` passed to the CV call |
| Feature selection | Top-k genes, highest variance or a univariate filter chosen on all samples, then CV of the classifier only | Selection inside the `Pipeline`; near-perfect CV on noise is the symptom |
| Label or proxy | A feature that is downstream of the outcome (post-diagnosis lab, treatment field, collection site tracking case/control); one feature dominating implausibly | Drop post-outcome variables; rerun without the dominant feature |
| Patient or replicate | Same patient, tumour, organoid or technical replicate in train and test; plain `KFold` | `routes/grouped.md` |
| Batch | Batch correlated with outcome; correction run across the split | Block the split by batch; never run unsupervised correction across it |
| Time | Random split of time-ordered data; future statistics used to standardise the past | `routes/site-or-time.md` |
| Duplicates | Near-identical samples, augmented copies, public-dataset overlap, homologous sequences | Deduplicate, or cluster then split, before CV |
| Test reuse | Repeated peeking to pick features, thresholds or "best epoch"; threshold set on the test fold | One locked test set; all tuning inside `routes/nested-cv.md` |

Other traps:

- `groups=` not threaded to the splitter, so it is ignored: pass it to `cross_validate` and `GridSearchCV.fit`.
- `cross_val_predict` used for the headline AUC: average per-fold scores, or pool per repeat as `scripts/cv_auc.py` does.
- A clean train/test split does not help if the leak happened before it.

Done when each row has been checked and the estimate was recomputed with the leaks fixed.
