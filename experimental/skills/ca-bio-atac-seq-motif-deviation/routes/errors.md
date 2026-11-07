# chromVAR problems

| Symptom | Cause | Fix |
|---------|-------|-----|
| `addGCBias` fails on missing `seqlengths` | Peak ranges lack chromosome sizes | `seqlengths(peaks) <- seqlengths(genome)` before the call |
| All z-scores near zero | Too few samples or near-identical samples | Needs 6 or more samples with biological variation; otherwise footprinting or differential accessibility |
| Variability all above 5, top motifs GC- or AT-rich | Under 5000 peaks or under 1500 reads per sample | Full-scale peakset; aggregate cells |
| Z-scores change between runs | Unseeded background sampling | Set and report `seed=` |
| `getBackgroundPeaks` slow | Large peakset | Keep defaults; do not reduce `niterations` below 30 |
| Every differential motif significant | No FDR applied, or wrong groups | Use `adj.P.Val`; check the `samples=` table |
| Signac `chromvar` assay all zero or NA | Motifs added before the peakset was final | Rerun `scripts/chromvar_signac.R` after the peaks are final |
| `FindAllMarkers` gives `avg_log2FC` on z-scores | Default fold-change function | Use `mean.fxn=rowMeans`, `fc.name='avg_diff'` (the script does) |
| Sign of the effect flipped | Reference level reversed | The first factor level is the reference; set `ref=` |
| ArchR `getMarkerFeatures` rejects NA | 0/0 z-scores in sparse cells | The script drops them and prints the count |
| Signac `RunChromVAR` not found | Removed in Signac 1.17 | Use `routes/signac.md` |
| `getMatrixSet` fails on the JASPAR2024 object | TFBSTools issue 39 | The scripts pass the SQLite handle |
| `plotVariability` warns about `aes_string()` | chromVAR with ggplot2 4.x | Harmless; the plot renders |

Done when the symptom's fix has been applied and the script rerun.
