# Common-Factor GWAS with Q_SNP
**Goal:** Identify SNPs that affect the latent factor and flag SNPs whose effect is heterogeneous across loadings.

**Approach:** Build SNP-by-trait effect-and-SE matrix via `sumstats()`, then fit the SNP-augmented model genome-wide via `commonfactorGWAS()`. Report both factor p-value and Q_SNP p-value per SNP.

`sumstats()` flips effect signs based on the reference panel A1/A2 to enforce consistent allele coding across input GWAS. Required input columns (case-sensitive) are `SNP, A1, A2, BETA/Z/OR, SE, P, N, MAF`; some are conditional on `se.logit` / `OLS`. Silent failures are almost always column-name mismatches (e.g. `EA`/`NEA` instead of `A1`/`A2` -> 0% SNPs retained) or a missing `N` column -> SNPs dropped. Always run `head(read.table(file, header=TRUE, nrow=2))` per input before the `sumstats()` call.

```r
# Prepare per-SNP betas and SEs across all input GWAS
ss <- sumstats(
    files = c('raw/trait1.txt', 'raw/trait2.txt', 'raw/trait3.txt'),
    ref = 'reference.1000G.maf.0.005.txt',
    trait.names = trait_names,
    se.logit = c(TRUE, TRUE, FALSE),     # TRUE if trait is logistic-scale (case-control); FALSE if continuous
    OLS = c(FALSE, FALSE, TRUE),
    linprob = c(FALSE, FALSE, FALSE),
    N = N,
    info.filter = 0.9,                   # standard imputation INFO threshold
    maf.filter = 0.01
)

saveRDS(ldsc_results, 'ldsc_results.rds'); saveRDS(ss, 'sumstats_snps.rds')
```

```bash
# commonfactorGWAS DWLS + ML cross-check + Q_SNP classification; last argument = cores
Rscript scripts/commonfactor_gwas_qsnp.R ldsc_results.rds sumstats_snps.rds cfgwas_qsnp.tsv 8
```

Output: the leading `sumstats()` columns (`SNP`, `CHR`, `BP`, `MAF`, `A1`, `A2` ...), then `i`, `lhs`, `op`, `rhs`, `est`, `se_c`, `Z_Estimate`, `Pval_Estimate` (factor effect), `Q` / `Q_df` / `Q_pval` (`Q_pval` IS the per-SNP Q_SNP test), `fail`, `warning`, and the script's `factor_sig`, `qsnp_sig`, `factor_only_dwls`, `Q_pval_ML`, `factor_only`. Cluster and `MPI=` settings: `references/runtime-and-python-tools.md`.

The "factor-only" subset (factor-significant AND Q_SNP non-significant) is the publication-grade set of common-factor SNPs.

**Cross-check Q_pval with `estimation = 'ML'`; DWLS Q_pval is provisional.** On two synthetic panels (flat SE; MAF/N-driven SE) built without `ldsc()`, `commonfactorGWAS(estimation = 'DWLS')` returned Q_pval 0.91-0.96 for both 5 planted heterogeneous SNPs and 5 planted factor SNPs (no separation), while `'ML'` on the identical input separated them (heterogeneous 6e-83 to 6e-11, factor 0.28-0.48; GenomicSEM 0.0.5 + lavaan 0.6.19, re-audit Inputs 2 and 8). Whether that is a DWLS weakness or an artifact of synthetic V without per-SNP N is unsettled. So: when V/N come from anything other than a genuine `ldsc()` + `sumstats()` run, treat DWLS Q_pval as provisional; and even on real data, do not call a SNP factor-only on a DWLS Q_pval alone.

The script runs that ML pass itself and sets `factor_only` only when DWLS and ML Q_pval both clear the threshold.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `usermodel()`/`commonfactorGWAS()`/`userGWAS()` error `object 'ReorderModel'` or `'ReorderModelnoSNP' not found` (both DWLS and ML) | lavaan >=0.7.0 requires `ordered=FALSE` for DWLS on continuous data; GenomicSEM 0.0.5 never supplies it (see Version Compatibility) | Pin lavaan to 0.6.19 (see Tool Installation); confirmed working under both estimators |
| `commonfactor()` reports "failed to converge" but `traceback()` shows an underlying `object '...Results' not found` error | Same lavaan-version incompatibility above, mislabeled by GenomicSEM's tryCatch as non-convergence | Not a model-specification problem -- do not re-specify the model; pin lavaan to 0.6.19 |
| `commonfactor()` complains "S not positive definite" | Genetic correlations near +/-1 among inputs | Drop redundant traits; verify rg < 0.95 pairwise |
| `ldsc()` fails with "category not found" | Wrong LD score column names (legacy format) | Use Python 3 LDSC fork; download `eur_w_ld_chr/` from alkesgroup |
| `lavaan` says "model not identified" | Too few traits for too many parameters | Need >= 3 traits per factor; constrain factor variance to 1 |
| MTAG `MaxFDR` not in log | Older MTAG version (< 1.0.7) | Update MTAG; MaxFDR reporting added late 2019 |
| `usermodel()` slow or fails | Complex syntax + many traits | Simplify model; estimate with `estimation = 'DWLS'`, not `'ML'`, when V is informative |
| Sumstats output has zero overlap with reference | Allele coding mismatch in `sumstats()` | Check `se.logit` and `OLS` settings per trait; align A1/A2 |
| GenomicSEM and TwoSampleMR give different rg | TwoSampleMR uses bivariate LDSC; GenomicSEM uses the same S | Match the underlying LDSC reference panel and weights |
