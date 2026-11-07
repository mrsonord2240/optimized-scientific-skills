# ArchR or paired multiome

```r
library(ArchR)
proj <- loadArchRProject('ArchR_out')
proj <- addCoAccessibility(proj, reducedDims='IterativeLSI', k=100, knnIteration=500, maxDist=250000)
co_acc <- getCoAccessibility(proj, corCutOff=0.5, returnLoops=TRUE)
```

The project needs a PeakMatrix; the loop result is the first element of the returned list in the source-tested version. See `references/other-methods.md` for return forms, LinkPeaks, SCENIC+ setup and published reference options. Paired RNA is needed for direct accessibility-expression linkage; a peak-pair score alone is not that linkage.

Done when the appropriate peak-pair or multiome output is produced and the method and cutoff are stated.
