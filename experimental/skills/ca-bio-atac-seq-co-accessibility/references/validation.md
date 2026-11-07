# Validation



```r
# Arc plot at ONE locus with a stringent cutoff: plotTracks(track) alone errors ("Unable to automatically
# determine plotting ranges"), and all strong arcs genome-wide draw an unreadable hairball.
library(Gviz); library(GenomicInteractions)
# Cicero Peak1/Peak2 are chr_start_end strings; convert to chr:start-end for GRanges()
to_gr <- function(x) GRanges(sub('_(\\d+)_(\\d+)$', ':\\1-\\2', x))
sel <- strong[strong$coaccess > 0.5, ]
gi <- GenomicInteractions(anchor1=to_gr(sel$Peak1), anchor2=to_gr(sel$Peak2),
                          counts=as.integer(sel$coaccess * 100))
track <- InteractionTrack(gi, name='co-accessibility > 0.5', chromosome='chr1', start=1e6, end=1.4e6)
plotTracks(list(GenomeAxisTrack(), track), chromosome='chr1', from=1e6, to=1.4e6)
```

For genome-browser visualization with ArchR: `plotPeak2GeneHeatmap()` shows the peak-gene linkage matrix; `plotBrowserTrack()` overlays connections on tracks.



For enhancer-to-gene linking with paired Hi-C/Micro-C, the canonical method is the ABC model (Fulco 2019, Nasser 2021), not Cicero. ABC computes ABC = (Activity_E * Contact_E,G) / sum_e(Activity_e * Contact_e,G); standardizes on combined ATAC + H3K27ac activity and Hi-C contact frequencies. ENCODE-rE2G (Gschwind et al 2023, bioRxiv) is the modern logistic-regression enhancer-gene link predictor.

See bio-atac-seq-enhancer-gene-linking for full ABC and ENCODE-rE2G coverage. Cicero is the ATAC-only fallback when no Hi-C is available.

## HiChIP H3K27ac as Orthogonal Anchor

| Decision | Action |
|----------|--------|
| Have Hi-C / Micro-C | Use ABC (bio-atac-seq-enhancer-gene-linking) primary; Cicero as ATAC-only sanity check |
| Have HiChIP H3K27ac | FitHiChIP loops (FDR < 0.05, count >= 5) primary; ABC + HiChIP intersection is high-confidence |
| Have ATAC + H3K27ac, no 3D | ABC with average HiC fallback (Fulco 2019); document degraded performance |
| Have only ATAC | Cicero (this skill); Hi-C concordance is partial |

Cicero is appropriate when no 3D data exists; do not use Cicero in lieu of ABC when Hi-C/Micro-C are available.

## Hi-C / Micro-C Concordance

Bands below are unsourced heuristics, not benchmarks; calibrate on your own Hi-C.

| Hi-C concordance | Action |
|-----------------|--------|
| > 50% of strong Cicero connections overlap Hi-C loops | High-confidence; Cicero captures real 3D structure |
| 30-50% | Common; some 3D contacts don't vary in accessibility |
| < 20% | Co-accessibility may not reflect contacts; lineage-specific contacts may be missing |

**Goal:** Quantify what fraction of strong Cicero connections are supported by Hi-C loop calls.

**Approach:** Import HiCCUPS loops as GenomicInteractions, build a parallel object from Cicero connections, then count anchor-anchor overlaps and report the percentage.

```r
# Compare Cicero against published Hi-C loops
library(GenomicInteractions)
hic_loops <- makeGenomicInteractionsFromFile('hiccups_loops.bedpe', type='bedpe',
                                             experiment_name='hiccups', description='HiCCUPS loops')
ci <- GenomicInteractions(anchor1=GRanges(sub('_(\\d+)_(\\d+)$', ':\\1-\\2', strong$Peak1)),
                          anchor2=GRanges(sub('_(\\d+)_(\\d+)$', ':\\1-\\2', strong$Peak2)))
overlap <- countOverlaps(ci, hic_loops) > 0   # anchor-anchor 'any' overlap; 'equal' is too stringent at loop bin resolution
cat(sprintf('Cicero connections overlapping HiCCUPS loops: %.1f%%\n',
            100 * mean(overlap)))
```

## Reconciliation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Cicero many weak connections; ArchR few strong | Different alpha or aggregation | Standardize parameters |
| LinkPeaks (Multiome) finds connections Cicero misses | LinkPeaks uses RNA expression as the anchor; Cicero is ATAC-only | Both valid; report intersection as high-confidence |
| Co-accessibility doesn't match Hi-C in heterochromatin | Heterochromatic contacts are constitutive; co-accessibility needs variation | Expected; co-accessibility complements Hi-C |
| SCENIC+ network has ENCODE-validated TFs but missing some | Motif database limited or RNA imputation missed | Expand motif database; integrate paired ChIP-seq if available |

**Operational rule:** Co-accessibility is a hypothesis generator. Validate with Hi-C, ChIP-seq, or experimental enhancer-promoter interaction (CRISPRi-FlowFISH).

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Cicero `make_cicero_cds` slow / crashes | k too high or cell count too large | Reduce k or subsample cells |
| All connections near zero | alpha set too high | Use `estimate_distance_parameter()` |
| Connection score > 1 reported | Bug in older Cicero versions | Update; check `as.numeric(coaccess)` for outliers |
| ArchR getCoAccessibility "TileMatrix" error | Need PeakMatrix not TileMatrix | `addPeakMatrix()` first |
| SCENIC+ install fails | Many heavy dependencies; pybedtools wheel build failed here | Use a Python <=3.11 env with conda-forge `pybedtools`; untested here; the upstream README documents no Docker image |
| Connection count varies wildly per run | Stochastic metacell aggregation | Set seed; or aggregate at higher k for stability |
| LinkPeaks all NaN | RNA expression has too many zeros | Re-filter cells with sufficient RNA |
| Peak names not matching | format mismatch (chr_start_end vs chr:start-end); the script accepts only `chr_start_end` | Standardize naming convention |
