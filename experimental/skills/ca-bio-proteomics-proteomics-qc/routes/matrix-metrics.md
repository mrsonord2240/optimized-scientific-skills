# Replicate correlation, sample swap, CV, missingness

```bash
python scripts/matrix_metrics.py proteinGroups.txt sample_annotation.csv --out matrix_metrics_results.tsv
```

Reads `LFQ intensity <sample>` by default (`--prefix` for another column set; correlation and CV need a quantified matrix, not the loading check). One TSV; the `table` column names each block: `replicate_correlation`, `cross_group_correlation`, `median_cv_linear`, `missingness_profile`, `completeness_filter` (`--min-valid-frac`, default 0.7 valid in at least one group).

- Summarize only `status == measured` rows. `not_measurable_n1` means unchecked, not clean. With no group of two or more the script stops.
- Correlation is on log2 (Pearson on raw intensity is a high-abundance artifact). Technical replicates r > 0.98; biological 0.90-0.98; below 0.8 a biological replicate is in doubt.
- A swapped pair still correlates 0.93-0.96 with its own group, so only `cross_group_correlation` shows a swap: `possible_swap = True` means another group matches better. It is a candidate; confirm against the sample sheet.
- CV is computed on the linear scale (SD/mean); on log2 data it comes out about 15-20x too small. Technical median CV < 10-20%, biological 20-40%. A lower CV is not automatically better. State transform, normalization and software with any CV.
- Missingness: present fraction rising with abundance is MNAR (left-censored); flat is MCAR. Match the imputer to it: kNN on MNAR invents mid-range values and kills present/absent calls; a left-shifted draw on MCAR fabricates false lows. Perseus-style down-shift is 1.8 SD, width 0.3. The imputation itself is `bio-proteomics-quantification`.
- Filter on completeness (valid in at least 50-70% of replicates in at least one condition) before imputing.

Done when `matrix_metrics_results.tsv` exists and each metric is stated for measured groups, with unmeasurable ones named.
