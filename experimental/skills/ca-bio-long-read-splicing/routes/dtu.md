# Differential transcript usage

```r
library(DRIMSeq); library(stageR)
```

Use transcript counts and a sample manifest with the pipeline in `references/dtu.md`. The source's FLAIR ID parsing assumes gene IDs contain no underscores; validate mapping for the data. Ambiguous and truncated reads leave assignment uncertainty.

Done when stage-wise gene and transcript results are reported.
