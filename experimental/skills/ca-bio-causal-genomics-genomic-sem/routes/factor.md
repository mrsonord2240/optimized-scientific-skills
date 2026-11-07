# Covariance and common factor

```bash
Rscript scripts/genomic_sem_commonfactor.R
```

First edit input paths, trait names, prevalence and sample sizes in the script. It also proceeds to GWAS, so prepare all inputs or stop after the factor stage. At least three indicators are needed for the source single-factor model. Inspect S and V, convergence, loadings, CFI and RMSEA; see `references/standard-workflow.md` and `references/failure-modes.md`. Use the compatible lavaan pairing.

Done when model fit, sampling covariance and limitations are reported.
