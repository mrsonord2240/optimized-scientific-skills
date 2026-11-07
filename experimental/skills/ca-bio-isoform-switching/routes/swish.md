# Uncertainty-aware transcript expression

```bash
Rscript scripts/swish_dte.R salmon_quant metadata.tsv dte_results.csv
```

Metadata TSV needs unique sample_id and condition. Quantifications need Gibbs or bootstrap replicates; the source uses Salmon --numGibbsSamples 20. This is DTE, not DTU. See `references/swish.md` for interpretation and uncertainty.

Done when transcript effects, q values and inferential uncertainty are saved.
