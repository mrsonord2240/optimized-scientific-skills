# DIA-NN run metrics (RT fit, peak width)

```bash
python scripts/diann_level1.py report.parquet --out diann_level1_results.tsv
```

Needs the main DIA-NN report (parquet or tsv) with `Run`, `RT`, `Predicted.RT`, `FWHM`, `Quantity.Quality`, `Global.Q.Value`; it refuses a matrix file. Per run, at `Global.Q.Value <= 0.01` (`--max-q`): `rt_fit_r2` (RT vs predicted RT), `median_fwhm` (minutes in DIA-NN 2.x), `median_quantity_quality` (reported only, no accepted cutoff), `fwhm_vs_median`, `flag`.

- `flag` is R^2 below 0.99 (`--min-r2`) or median FWHM more than 1.25x the across-run median (`--max-fwhm-ratio`; the accepted band is a 20-30% rise over baseline).
- The across-run median stands in for a rolling baseline, and with three runs one broad run shifts it. The script warns under 4 runs. Trend over the queue (Levey-Jennings) before acting on one flag.
- For protein counts, require both `Global.Q.Value` and `Global.PG.Q.Value` at or below 0.01: precursor q is not protein q. Internals are `bio-proteomics-dia-analysis`.

Done when `diann_level1_results.tsv` exists and flagged runs are named with the value that fired.
