# Independent-group feature tests

```bash
python scripts/differential_csv.py peaks.csv samples.csv case control results.csv
```

The peak CSV has a `feature` column and sample columns; metadata CSV has unique `sample,group` columns. Every peak sample must be present in metadata. The wrapper follows the source's log2 transform, Welch tests with at least three observed positive values per group, difference of log means and BH adjustment across valid p-values. Zeros and missing intensities are excluded; negatives are rejected. It does not support pairing, covariates, missingness modeling or effect confidence intervals; use a suitable design-specific model for those. Output may be overwritten. Check detection rates by group and dependence assumptions.

Done when the experimental unit, test family, effect direction, missingness and uncertainty are reported.
