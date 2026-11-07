# Iterative overlap, 501 bp fixed width

```bash
bash scripts/iterative_overlap.sh 'peaks/per_sample/*.narrowPeak' hg38.chrom.sizes hg38-blacklist.v2.bed.gz 250 consensus_out
```

Arguments in order: quoted narrowPeak glob, chrom.sizes, blacklist (BED or .gz), half-width, output directory. Review the paths: defaults are hg38 names.

The script validates every row (ten tab-delimited columns, usable summit offset; `-1` is rejected with file and line), pools the peaks, re-centers each on its summit, extends to `2*half_width+1` bp, sorts by signalValue descending, keeps each peak that overlaps no kept peak, removes blacklist overlaps, and writes a SAF.

- Output: `consensus_out/consensus_final.bed` (BED), `consensus_out/consensus.saf` (featureCounts), plus intermediates.
- Do not replace the greedy step with `bedtools cluster` plus one peak per cluster: clusters chain transitively and discard peaks the greedy rule keeps.
- A smaller half-width (for example 100) keeps more peaks when many summits lie within 500 bp.
- Single-cell ATAC: pool MACS3 calls per cluster and pass them all as inputs.

Done when `consensus_final.bed` exists; next `routes/count.md`. R only: `routes/iterative-overlap-r.md`.
