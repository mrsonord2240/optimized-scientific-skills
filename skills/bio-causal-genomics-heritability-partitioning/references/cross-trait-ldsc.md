# Cross-trait LDSC genetic correlation (reference for heritability-partitioning SKILL.md)

## Cross-Trait LDSC for Genetic Correlation

**Goal:** Estimate genetic correlation `rg` between two traits with calibrated handling of sample overlap.

**Approach:** Munge both sumstats with identical SNP list -> run --rg with two munged files; the bivariate intercept absorbs sample overlap and the rg estimate remains unbiased.

```bash
ldsc.py \
    --rg trait1.sumstats.gz,trait2.sumstats.gz \
    --ref-ld-chr eur_w_ld_chr/ \
    --w-ld-chr eur_w_ld_chr/ \
    --out rg
# Output: rg, se, p, gcov_int (cross-trait intercept), h2_obs, h2_int per trait
```

The cross-trait intercept (`gcov_int`) is the LDSC analog of sample-overlap z-score correlation; non-zero indicates sample overlap or cryptic shared structure. LDSC rg is unbiased even with sample overlap because the bivariate intercept absorbs it. HDL is more precise but requires non-overlapping samples.

Non-zero `gcov_int` under known sample overlap is the **correct** behavior, NOT pathology. The rg estimate remains unbiased; the intercept is the overlap absorber, doing its job. Pre-empt the reviewer comment "gcov_int = 0.05 with a shared cohort is expected, not confounding evidence" by reporting `gcov_int` alongside rg and explaining the absorber role.
