# Iterative overlap in R

```bash
Rscript scripts/iterative_overlap.R hg38.chrom.sizes hg38-blacklist.v2.bed.gz consensus.bed 250 rep1.narrowPeak rep2.narrowPeak rep3.narrowPeak
```

Same algorithm as `routes/iterative-overlap.md` with GenomicRanges: validated narrowPeak rows, summit re-centering to `2*half_width+1` bp, greedy non-overlap by descending signalValue, blacklist removal (pass `NONE` to skip). It writes `consensus.bed` only; use `routes/count.md` for SAF conversion.

- `!duplicated(findOverlaps(...))` or `reduce()` or `disjoin()` do not implement the greedy rule. Keep the loop.
- The loop is quadratic. Above about 500k pooled peaks use `routes/iterative-overlap.md`.
- `trim()` clamps to chromosome ends only because the script sets `seqlengths` from chrom.sizes.

Done when `consensus.bed` exists; next `routes/count.md`.
