# HyPrColoc

## HyPrColoc Cluster-Based Coloc (Many Traits)

HyPrColoc (Foley 2021) extends single-causal coloc to many traits by clustering traits that share a causal variant. Output: per-cluster posterior + per-trait cluster assignment.

```r
library(hyprcoloc)
# Inputs: SNPs-by-traits matrices of betas and standard errors
# Rows = SNPs (must be shared across all traits); Columns = traits
res <- hyprcoloc(effect.est=betas, effect.se=ses,
                  trait.names=colnames(betas), snp.id=rownames(betas),
                  reg.thresh=0.7,     # regional probability of coloc threshold
                  align.thresh=0.7)   # alignment threshold for traits within a cluster
res$results  # cluster assignment per trait + posterior
```

HyPrColoc inherits the single-causal-per-cluster assumption from coloc.abf; clusters can fragment if the true biology is allelic heterogeneity.
