# ArchR and SnapATAC2

```r
library(ArchR)
addArchRGenome('hg38')
```

For complete fragment import, QC, doublet filtering, LSI and peak-matrix code use `references/ecosystem-workflows.md`. SnapATAC2 uses `pp.import_fragments` in the source-tested version; `pp.knn` precedes Leiden. Select genome and annotations matching the data, use TileMatrix for embedding and PeakMatrix for downstream peak analysis, and inspect empty Arrow results.

Done when the requested project, clustering and relevant peak matrix are produced and QC recorded.
