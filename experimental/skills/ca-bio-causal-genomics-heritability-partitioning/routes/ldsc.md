# Total and partitioned LDSC

```bash
bash scripts/ldsc_partitioned_h2.sh gwas.tsv trait 0.08 0.10
```

The optional prevalences are sample case fraction and population prevalence, respectively; omit both for quantitative traits. Supply the reference bundles via LDSC_REF and match ancestry. The wrapper runs total, baseline-LD and cell-type analysis; source CBIIT version needs the documented h2-cts patch in `references/compatibility.md`. Adjust the workflow if cell-type analysis is not requested. Correct annotation tests across the declared family, report enrichment and tau, and check model assumptions using `references/interpretation.md`.

Done when requested outputs, scale, diagnostics and uncertainty are reported.
