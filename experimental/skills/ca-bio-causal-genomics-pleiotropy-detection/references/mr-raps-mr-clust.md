# MR-RAPS and MR-Clust

## MR-RAPS Loss Function and Overdispersion

**Trigger:** Weak instruments (mean F < 20) and/or suspected UHP requiring outlier-resistant estimation.

- `over.dispersion = TRUE` always for MR (horizontal-pleiotropy variance is real, not noise; turning this off underestimates SE)
- `loss.function = 'huber'` (default; outlier-resistant; suited to mild to moderate UHP)
- `loss.function = 'tukey'` (more aggressive; downweights extreme outliers more; choose when many obvious outliers suspected)
- `loss.function = 'l2'` (non-robust; equivalent to weighted least squares; do not use when UHP suspected)

These are not top-level arguments to `TwoSampleMR::mr_raps()` -- its signature is
`mr_raps(b_exp, b_out, se_exp, se_out, parameters = default_parameters())`, so pass them nested:

```r
TwoSampleMR::mr_raps(b_exp = dat$beta.exposure, b_out = dat$beta.outcome,
                      se_exp = dat$se.exposure, se_out = dat$se.outcome,
                      parameters = list(over.dispersion = TRUE, loss.function = 'huber', shrinkage = FALSE))
```

Calling with bare `over.dispersion = TRUE, loss.function = 'huber'` throws `unused arguments`.

Tukey is preferable when leave-one-out reveals 2+ SNPs single-handedly shifting the IVW estimate by > 1 SE.

## MR-Clust for Mechanism Heterogeneity

**Goal:** When a single causal estimate is misleading because instruments operate through multiple causal mechanisms (e.g. LDL on CHD via multiple lipoprotein subfractions), identify clusters of instruments with similar per-SNP Wald ratios.

```r
library(mrclust)
ratio_hat <- dat$beta.outcome / dat$beta.exposure
ratio_se <- abs(dat$se.outcome / dat$beta.exposure)
res_mc <- mr_clust_em(theta=ratio_hat, theta_se=ratio_se,
                     bx=dat$beta.exposure, by=dat$beta.outcome,
                     bxse=dat$se.exposure, byse=dat$se.outcome,
                     obs_names=dat$SNP)
per_cluster <- res_mc$results$best
```

Clusters with cluster_class = 'null' are pleiotropy-only instruments. Per-cluster IVW estimates may differ substantially; biological annotation of the SNPs in each cluster (pathway, target gene) is the interpretation step. Pure statistical clustering without a biological story is weak evidence.
