# Drift, normalization and imputation

```bash
Rscript examples/normalize_data.R
```

This is a synthetic base-R demonstration that writes and deletes a temporary result; it is not a real-data pipeline. Production QCRSC, PQN, batch correction and QRILC/missForest paths require input-specific adaptation and installed-package verification. Validate QC coverage, exclude conditioning injections and inspect order/batch confounding before correction. Hold out QCs and inspect dilution linearity and biological variance; tighter fitted-QC clustering alone is insufficient. Filter detection before imputation; QRILC requires a positive log-scale round-trip and explicit missingness assumptions. Compare plausible imputation choices. See [assumptions](../references/assumptions.md).

Done when preprocessing choices, exclusions, held-out validation and sensitivity to assumptions are reported.
