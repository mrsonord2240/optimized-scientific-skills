# Raw signal and contaminants (run before normalizing)

```bash
python scripts/raw_qc.py proteinGroups.txt sample_annotation.csv --out raw_qc_results.tsv
```

`sample_annotation.csv` has a `sample` column (the text after `Intensity `) and a `condition` column; other names: `--sample-col`, `--group-col`. The script removes `Potential contaminant`, `Reverse` and `Only identified by site` rows, then writes one row per sample: `n_quantified`, `total_signal`, `contaminant_pct` (share of summed raw intensity, measured before the rows are removed), `fold_total_vs_group`, `ids_vs_group`, `flag`, `loading_rule`. It refuses an `LFQ` or `Normalised` prefix.

- Loading rule: total signal at or below 0.5x the group median, or ID count below 0.8x it, is a loading or injection failure to exclude, not rescale. Judge on total signal and ID count, not the boxplot median; left-censoring hides a low load there. Override with `--min-fold-total` and `--min-ids-frac` only when the request gives its own.
- `loading_rule = fallback_all_samples` means a group had one sample and the within-group rule could not fire. Report that, and that a loading difference tracking condition is not separable from biology in that design.
- Contaminant fraction has no universal cutoff: compare against the lab baseline for that sample type and between groups. Keratin and trypsin dominate low-input samples.
- A flagged sample goes in the exclusion log (`routes/writeup.md`). Do not drop it silently.
- DIA-NN input: use `Precursor.Quantity` per run, not the `PG.MaxLFQ` matrix.

Done when `raw_qc_results.tsv` exists, flagged samples are listed with value and threshold, and the contaminant fraction is stated per group.
