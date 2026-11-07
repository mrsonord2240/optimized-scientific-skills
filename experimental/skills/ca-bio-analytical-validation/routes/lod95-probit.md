# LoD95 and LoB from a dilution series

```bash
python scripts/lod95_probit.py dilution.csv --blanks blanks.csv --out lod95_results.tsv
```

`dilution.csv`: either one row per replicate (`vaf`, `detected` as 0 or 1) or one row per level (`vaf`, `n_replicates`, `n_detected`); other column names: `--vaf-col`, `--detected-col`. `blanks.csv` (optional) has a `signal` column of blank-replicate readings on the VAF scale. The script fits a probit GLM of detection on log10(VAF) and writes a metric/value table: `lod95_vaf`, its confidence interval (`--confidence`, default 0.95), slope, convergence, deviance, the empirical detection rate per level, and `lob_gaussian` (mean + 1.645 SD of the blanks, at least two needed).

The script refuses, rather than return a precise-looking number: fewer than three VAF levels, fewer than `--min-replicates` (default 3) at any level, empirical rates that do not bracket 0.95, perfect separation, non-convergence, a non-positive slope, an LoD95 outside (0, 1], or a VAF above 1. Report the refusal as the result and ask for replicates near the transition; do not loosen the checks.

- Report LoD95 with its interval and slope, never the point alone, and condition it on input mass, consensus depth and replicate count: "LoD95 0.1% VAF at 30 ng, 2x duplex, 95% of 20 replicates".
- LoD is the yes/no call; LoQ (for example CV at or below 20%) is separate and is never below LoD. Below LoQ report detected or not detected, not a VAF value.
- The Gaussian LoB is a shortcut; a full EP17 study may need finite-sample or nonparametric treatment. In plasma, clonal haematopoiesis variants dominate the non-tumour blank signal.
- The fit takes pooled replicate binary outcomes. A teaching-model dilution series (simulated from a chosen true LoD) is not validation evidence.

Done when `lod95_results.tsv` exists, or the refusal is reported with what data would fix it.
