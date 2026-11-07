# Uncertainty-aware transcript expression

```r
library(fishpond); library(tximeta)
```

Salmon quantifications need Gibbs or bootstrap replicates; use `--numGibbsSamples 20` at quantification. Construct the metadata and run `references/swish.md`; condition must be a factor. This is DTE, not a DTU test.

Done when transcript effects, q values and inferential uncertainty are reported.
