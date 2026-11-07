# Two-group comparison

```bash
Rscript scripts/consensus_csv.R counts.csv metadata.csv consensus.csv
```

Counts CSV has taxon IDs in the first column and samples in columns; metadata CSV has `sample,Group` with exactly two groups and one row per sample. The wrapper aligns IDs, requires nonnegative integer counts, filters taxa present in at least 10% of samples, and follows the source ALDEx2 settings (128 draws, all-feature denominator, expected BH < 0.05 and absolute effect > 1) plus LinDA at alpha 0.05. It writes the union with method membership and overwrites the output. This is a cross-sectional workflow and does not handle repeated subjects or covariates. Inspect full per-method results and sensitivity to filtering and reference assumptions; the wrapper only saves the hit union.

Done when sample alignment, design, prevalence threshold, direction, per-method evidence and disagreement are reported.
