# Stage-wise differential transcript usage

```r
library(tximport); library(DRIMSeq); library(DEXSeq); library(stageR)
```

Supply named Salmon quantification files, metadata and transcript-to-gene mapping, then run `references/manual-dtu.md`. It filters with DRIMSeq, tests with DEXSeq, and uses stageR for gene and transcript decisions. For large groups consider satuRn as described in `references/selection.md`.

Done when gene-level and transcript-level adjusted results are available.
