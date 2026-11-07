# Is a split real? Significance and valid marker inference

Marker tests (`rank_genes_groups`, `FindAllMarkers`) after clustering are conditioned on labels chosen to maximize separation. Under one homogeneous population their type-I error approaches 1 as resolution rises, and multiple-testing correction does nothing. Use them to rank and label only.

For a defensible claim that clusters are distinct:

1. A cluster significance test. `scSHC` takes an HVG expression matrix, genes in rows and cells in columns; a failed split is over-clustering, so merge. CHOIR is the alternative. Both are optional installs:

```r
if (!requireNamespace('scSHC', quietly = TRUE)) {
  stop("Install scSHC from github.com/igrabski/sc-SHC before making a cluster-significance claim")
}
# expr_hvg: normalized or transformed HVG matrix, genes x cells; choose preprocessing options per the package docs
significant_clusters <- scSHC::scSHC(expr_hvg)
```

2. A double-dipping-robust DE method for the markers: ClusterDE (synthetic null), or count splitting where the noise model holds. Poisson thinning of overdispersed counts silently reinstates the bias.
3. Replication in an independent sample.

If scSHC is not installed, say so and report the split as unvalidated. Do not substitute marker p-values.

For cross-condition DE use the pseudobulk route of the DESeq2 Skill; for marker testing see markers-annotation.

Done when the claim states which test supported it, or states that none was run.
