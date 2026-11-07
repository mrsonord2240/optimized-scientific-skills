# ArchR scATAC

```bash
Rscript scripts/chromvar_archr.R proj_dir markers.rds group=Clusters motifset=cisbp
```

For a project without a peak set add `peakset=build macs2=/path/to/macs2`. The script runs `addPeakMatrix`, `addMotifAnnotations`, `addBgdPeaks`, `addDeviationsMatrix`, drops cells with NA z-scores (printing how many) and runs `getMarkerFeatures` on `MotifMatrix` with `useSeqnames='z'`.

- `motifset=cisbp` is ArchR's default (870 motifs in the ArchR 1.0.3 test); `JASPAR2020` is the alternative. Choose one and report it.
- Use the PeakMatrix, not tiles: tile deviations add intergenic noise and are for embedding.
- `getMarkerFeatures` refuses NA z-scores. In the test slice they came from cells with 26-203 peak reads whose rare motifs had no matched reads. Count them first, report how many cells were dropped, and do not assume the same on a full-depth dataset.
- Under 500 cells per cluster z-scores may be unreliable.

Done when `markers.rds` exists and the answer reports the number of cells dropped and the motif set.
