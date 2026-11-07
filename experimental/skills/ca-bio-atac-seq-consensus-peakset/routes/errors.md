# Consensus peakset problems

| Symptom | Cause | Fix |
|---------|-------|-----|
| Differential peaks track peak width | Variable-width consensus | Fixed width: `routes/iterative-overlap.md` or `routes/diffbind.md` |
| Near-zero significant differential peaks | Intersection across all samples dropped condition-specific peaks | `routes/per-condition.md` |
| Script stops with "summit offset is -1" or a column-count error | narrowPeak with no summit, or not 10 columns | Use peaks called with summits; broadPeak cannot be summit-centered |
| Negative coordinates | Summit near chromosome start | `bedtools slop -g chrom.sizes -b 0`, or `awk '$2 >= 0'` |
| Output much smaller than expected | Many summits within 500 bp, strict non-overlap | Smaller half-width (for example 100) |
| `multiinter` output empty | Unsorted inputs or mismatched chromosome naming | Sort every input; make names consistent |
| featureCounts non-unique GeneID | Duplicate peaks in the SAF | Non-redundant set; `chr_start_end` IDs |
| Peak counts differ across builds | Liftover loss, build-specific blacklists | `routes/liftover.md`; match the blacklist to the build |

Done when the symptom's fix has been applied and the route rerun.
