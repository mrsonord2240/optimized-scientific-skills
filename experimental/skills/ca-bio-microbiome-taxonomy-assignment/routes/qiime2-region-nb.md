# QIIME2 region-matched naive Bayes (classify-sklearn)

```bash
F_PRIMER=GTGYCAGCMGCCGCGGTAA R_PRIMER=GGACTACNVGGGTWTCTAAT \
  bash scripts/qiime2_region_nb.sh silva-138-99-seqs.qza silva-138-99-tax.qza rep-seqs.qza taxonomy_out
```

Primers shown are 515F/806R (V4); use yours, and give the reverse primer 5' to 3', not reverse-complemented. Three steps, in this order:

1. `extract-reads` cuts the reference to the primer-bounded region (`--p-min-length 50 --p-max-length 0`).
2. `fit-classifier-naive-bayes` trains on the extracted reads.
3. `classify-sklearn --p-confidence 0.7 --p-read-orientation auto --p-n-jobs 1` writes `taxonomy_out/taxonomy.qza`.

`CONFIDENCE` and `N_JOBS` env override 0.7 and 1.

- Downloading a region-matched pre-trained `.qza` is the shortcut, but it must be built for this exact QIIME2 release (`routes/errors.md`); otherwise train locally.
- Confidence 0.7: below it the lineage is truncated to the deepest rank that clears it. 0 computes but never truncates; `disable` skips. State the value used.
- Training memory: a full reference (SILVA or GTDB, 400K+ sequences) needs tens of GB of RAM and can be OOM-killed; subsample the extracted reference (about 60,000 sequences) or use a smaller one.
- `--p-n-jobs` above 1 multiplies memory, one classifier copy per job; a full-length SILVA classifier is multi-GB. Lower `--p-n-jobs` or `--p-reads-per-batch` if OOM-killed.
- Species-aware option for a known habitat (gut, soil): a weighted (clawback) classifier with a habitat abundance prior raises species accuracy.

Done when `taxonomy.qza` exists and the database release, primer region and confidence are recorded.
