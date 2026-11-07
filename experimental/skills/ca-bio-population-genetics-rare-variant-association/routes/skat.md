# Direct SKAT

```r
obj <- SKAT::SKAT_Null_Model(phenotype ~ age + sex + PC1 + PC2, out_type = 'D', data = covar_df)
skato <- SKAT::SKAT(Z, obj, method = 'SKATO', weights.beta = c(1, 25))
```

This is an adaptation snippet, not a complete file-based runner. Supply a prepared covariate/phenotype data frame and sample-by-variant genotype matrix in exactly the same sample order; choose binary (`D`) or continuous (`C`) null as appropriate. The source also uses `r.corr=1` for burden and the default SKAT test for the variance-component endpoint. For genome-wide SSD, compare the `.fam` IID order to the covariate IID order before fitting: matching is by position. Relatedness and severe case-control imbalance may require a suitable calibrated null beyond this small-cohort example. Confirm installed SKAT interfaces.

**Done when:** p-values, genotype/sample alignment, weights, masks, null covariates and test family are documented; build a file-loading/export wrapper appropriate to the data before execution.
