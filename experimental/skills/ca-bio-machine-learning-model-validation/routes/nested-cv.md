# Tuned pipeline: nested cross-validation

```bash
python scripts/nested_cv.py data.csv --label label --drop sample_id --ks 10,50,200 --Cs 0.01,0.1,1
```

Pipeline: scaling, `SelectKBest`, logistic regression. The inner loop chooses `k` and `C`; the outer loop grades the winner once on folds it never saw. Add `--group patient_id` when rows share a patient (both loops then split by patient). Defaults: 5 outer x 5 inner folds, `--repeats 1`.

- The nested AUC is the estimate. Also state the flat score the script prints (best inner score on all data), the number of configurations searched, and n. The gap is the optimism.
- Nested CV is needed for any choice made by looking at results, including "tried three options, kept the best".
- Pass your own grid by editing `--ks` and `--Cs`; a larger grid means larger optimism.
- Never report the flat score as the performance.

Done when you have reported the nested AUC (mean, SD over outer folds) and the flat score beside it.
