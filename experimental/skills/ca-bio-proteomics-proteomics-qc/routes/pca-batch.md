# PCA and batch

```bash
python scripts/pca_batch.py proteinGroups.txt sample_annotation.csv --batch-col batch --out pca_batch_results.tsv
```

Run on the normalized survivors of `routes/raw-signal.md`. Default input is `LFQ intensity <sample>` (`--prefix` to change). The TSV has one row per sample (PC coordinates, annotation columns) plus `__explained_variance__` and `__batch_test_PCn__` rows (`p`, `status`). It uses complete proteins only and a seeded exact solver, so p-values repeat.

- It refuses a design where every batch level holds a single condition: batch and condition are then the same variable and nothing can be corrected or tested. Report it as non-identifiable and stop.
- `status = not_testable` (a batch level with one sample) is never evidence of no batch effect.
- If a top PC tracks batch, keep batch in the design for the differential test. Use `limma::removeBatchEffect` or ComBat only on the matrix used for PCA and plots; do not test on a corrected matrix and also model batch.
- A sample isolated from its group is a removal or re-run candidate, but a high-missing sample is judged by `raw_qc.py`, not by PCA.
- Stop and ask before excluding samples, when n < 5 per group (PCA is unstable), or when no un-normalized column exists for the loading check. Re-run the downstream check with and without borderline samples.

Done when the PC variance, the batch test status per PC and any outlier candidates are stated, and confounding is addressed.
