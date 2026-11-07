# Signac scATAC

```bash
Rscript scripts/chromvar_signac.R in.rds out.rds markers.csv assay=ATAC
```

`in.rds`: a Seurat object with a ChromatinAssay on the final peakset and cluster identities set as the active identity. The script runs `AddMotifs` (JASPAR2024 vertebrates CORE), drops empty peaks and peaks with undefined GC, samples background peaks under a fixed seed, runs chromVAR directly on the peak counts, stores z-scores as a `chromvar` assay and runs `FindAllMarkers(only.pos=TRUE, mean.fxn=rowMeans, fc.name='avg_diff')`. It writes the marker table with `motif_name` added.

- `Signac::RunChromVAR()` is gone in Signac 1.17 (present in 1.16.0 and earlier). The script does not use it.
- `mean.fxn=rowMeans` is required: the default fold-change is for log counts and is meaningless on z-scores. Marker row names are JASPAR IDs.
- NA z-scores: the script prints the count. NAs come from motifs added before the peakset was final, or sparse cells. Drop them and report how many; fix the cause rather than ignoring them.
- Under 500 cells per cluster z-scores may be unreliable; aggregate cells or enlarge the peakset.

Done when `markers.csv` and `out.rds` exist and the answer reports NA count, seed and the motif set.
