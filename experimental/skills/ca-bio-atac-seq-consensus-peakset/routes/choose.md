# Choosing a strategy

| Goal | Strategy | Route |
|------|----------|-------|
| Modern ATAC publication, ML features, cross-study comparison | Iterative overlap, 501 bp fixed width | `routes/iterative-overlap.md` |
| Standard 2-3 replicate differential with DiffBind | DiffBind summit re-centering, `minOverlap=2` | `routes/diffbind.md` |
| Multiple conditions with a strong biology shift | Consensus per condition, then union | `routes/per-condition.md` |
| ENCODE-compliant reproducibility | Union of IDR-passed peaks from true replicate pairs | `routes/union.md` |
| Single-cell ATAC pseudobulk | MACS3 per cluster, then iterative overlap across clusters | `routes/iterative-overlap.md` |
| Cross-study, other build | Iterative overlap on one build after liftover | `routes/liftover.md` |
| Quick look only | `bedtools merge` (width-biased; never for statistics) | `routes/union.md` |

Why not the others: a naive union of variable-width peaks makes differential results track peak width; requiring a peak in all replicates drops condition-specific peaks; a global majority rule drops a peak present in all replicates of the smaller condition (2 of 5 is under 50%); extending from the peak midpoint instead of the summit shifts off the biology.

Two special cases:

- ATAC plus broad histone marks (H3K27ac, H3K4me1) in one study: a 501 bp window loses most of a 2-5 kb domain. Build 501 bp for ATAC and a broadPeak union for the histone mark, and run differential on each separately.
- ENCODE-rE2G predicted regulatory elements can stand in as a consensus for ML features, but only for ENCODE cell types and they may not align to your dataset's peaks.

Done when one strategy is chosen and its route is open.
