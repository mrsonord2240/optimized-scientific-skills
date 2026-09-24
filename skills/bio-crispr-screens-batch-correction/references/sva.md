## SVA (Surrogate Variable Analysis)

**Goal:** Estimate unknown latent factors that may confound the screen.

**Approach:** SVA computes surrogate variables that capture variance not explained by known biological factors; these can then be added to the MAGeCK MLE design matrix as covariates.

```r
library(sva)
# counts_df: rows = sgRNAs, columns = samples
mod <- model.matrix(~ condition, data = metadata)
mod0 <- model.matrix(~ 1, data = metadata)
sv_obj <- sva(as.matrix(counts_df), mod, mod0)
n_sv <- sv_obj$n.sv  # number of surrogate variables
# Add to design matrix for MAGeCK MLE
design_mat <- cbind(mod, sv_obj$sv)
```

**Use case:** When the screen has clear biological signal (e.g. essentiality recovery passes) but small effect sizes are hidden by noise; SVA-discovered latent factors as covariates can recover them.
