# New site, or later samples

```bash
python scripts/split_auc.py data.csv --label label --drop sample_id --site site
python scripts/split_auc.py data.csv --label label --drop sample_id --time collection_date
```

`--site`: leave-one-site-out; each site is the test set once and the AUC is reported per site. `--time`: rows are sorted by the column and split forward in time (`TimeSeriesSplit`, `--splits 5`); training data always precedes the test block. Pipeline as in `routes/fixed-pipeline.md` (`--k`, `--C`), refit per split.

- Use exactly one of `--site` and `--time`.
- Never shuffle time-ordered data or fit statistics on later periods.
- A site or split with one class has no AUC; the script skips it and says so.
- Report the per-site (or per-split) AUCs, their mean, and their range. A wide range is the finding.
- If batch tracks the outcome, do not run batch correction across the split.

Done when you have reported per-site or per-split AUCs and their mean.
