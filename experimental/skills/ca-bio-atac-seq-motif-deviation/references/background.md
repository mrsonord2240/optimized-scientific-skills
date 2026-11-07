# Background (not needed to run an analysis)

Tools:

| Tool | Background | Notes |
|------|------------|-------|
| chromVAR | Matched on GC and accessibility, 50 peaks per match | Bulk and single cell; sparse-aware |
| Signac + chromVAR | Same | `RunChromVAR` only in Signac 1.16 and earlier |
| ArchR `addDeviationsMatrix` | `getBgdPeaks`, GC and log accessibility | ArchR-specific format; faster on large data |
| SCENIC+ | Multi-modal | Needs paired RNA; chromVAR alone is not enough |

Background matching: each foreground peak is matched on GC and total accessibility using bin size `bs` (default 50), sampling `niterations` (default 50) replacements; the variance across these samples is the null. Below 30 iterations the noise inflates variability; above 100 gains diminish. For very small peaksets (under 2000) lower `bs` to avoid empty bins. For non-canonical genomes consider rebuilding bins with `quantile()`.

Why correction compresses broadly accessible cell types: the correction normalizes total accessibility, so high-background types get compressed z-scores and PCA on z-scores separates cell types less cleanly than raw counts. Hence per-cluster runs when accessibility differs more than 5x.

Heuristics inherited from the source, not tested or cited: 5000 peaks, 500 cells per cluster, the 5x accessibility rule, the 3-5 versus 6 or more sample rules, and chromVAR plus scBasset intersection. Only the 6-sample bulk run and the two single-cell chains were executed; time-course splines, custom PFMs and scBasset were not. Example counts: GM12878 versus K562 gave top |z| about 7 and 742 of 879 motifs passing `|logFC| >= 0.5`; JASPAR2024 CORE vertebrates has 879 motifs (1,912 with `all_versions=TRUE`). Unseeded reruns changed z-scores by up to 1.25 and the significant set by about 12%.

Sequence-model alternatives: scBasset (sequence CNN with per-cell projection; higher cluster discrimination, newer, needs about 100 cells per cluster) and Enformer-derived TF activity (long context, cell-type-specific pretrained models).

The bulk script in the source hard-coded a 6-sample condition vector and the hg38 genome; `scripts/chromvar_bulk.R` takes them as arguments. The Signac and ArchR code in the source reference file became `scripts/chromvar_signac.R` and `scripts/chromvar_archr.R`. The source tested those on a 10x PBMC 5k chr1:1-30 Mb slice. In the ArchR slice the NA-z cells had 26-203 peak reads (median cell 917).
