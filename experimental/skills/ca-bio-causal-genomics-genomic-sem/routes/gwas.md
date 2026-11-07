# Factor GWAS and Q_SNP

```bash
Rscript scripts/commonfactor_gwas_qsnp.R ldsc_results.rds sumstats_snps.rds cfgwas_qsnp.tsv 8
```

Inputs are the full LDSC covariance object and harmonized sumstats output. The wrapper runs DWLS and ML and calls factor-only only when both Q tests clear 0.05 divided by factor-significant SNP count. DWLS Q alone is provisional; source comparisons found unresolved estimator differences. Inspect failures and convergence; see `references/gwas-and-errors.md`.

Done when factor effects and Q_SNP results, cross-check and assumptions are reported.
