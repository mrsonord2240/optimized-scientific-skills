# Fixed pipeline, one sample per patient

```bash
python scripts/cv_auc.py data.csv --label label --drop sample_id --k 50
```

`data.csv`: one row per sample, `label` column (binary), `--drop` for ID columns, every other numeric column is a feature. The pipeline is univariate selection of `--k` features, scaling, then logistic regression (`--C`), all refit inside each training fold. Defaults: 5 folds x 10 repeats.

- The estimate is the mean AUC over repeats; its SD is the uncertainty. Give both and stop.
- Do not shop between 5-fold, 10-fold and leave-one-out for the number. Leave-one-out cannot give an AUC.
- Per-fold AUCs on a dozen samples are noisy; the script pools out-of-fold scores within a repeat and computes one AUC per repeat.
- The script warns if a `--drop` column has repeated values. If rows share a patient or donor, switch to `routes/grouped.md`.
- If `--k` or `--C` were picked by looking at results, the number is optimistic: use `routes/nested-cv.md`.
- To validate your own model, import `cv_auc` and call `repeated_cv_auc(pipe, X, y, groups)` with a `Pipeline` that holds every data-dependent step.

Done when you have reported mean AUC and SD from the script's output line.
