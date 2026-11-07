# Discriminant model validation

```bash
Rscript examples/metabolomics_stats.R
```

The source runs three synthetic ropls demonstrations with 1,000 permutations by default and may take substantial time; it is not a user-data CLI. Adapt the guarded fit to a real matrix and design, and record whether OPLS-DA or its PLS-DA fallback was fitted. Empty model summaries must not be read as results. A supervised score plot, R2 or VIP > 1 alone is not evidence of predictive performance; perform full-pipeline permutation and nested validation with preprocessing and selection inside folds, and seek independent validation. Compare Pareto and unit-variance sensitivity. No generic covariate, paired, mixed-model or leakage-safe external validation workflow is provided here.

Done when model choices, Q2, permutation design, validation and uncertainty are reported without treating a score plot as proof.
