# ENCODE-style TSS enrichment

```bash
bamCoverage -b filtered.bam -o sample.bw -bs 1 --extendReads --normalizeUsing None
python scripts/encode_tss_enrichment.py sample.bw tss.bed6 --flank 2000 --output score.json
```

The score is mean signal in the 100 bp window at the TSS over mean signal in the outer 100 bp at each end of the +/-2000 bp window, on the average TSS profile. The script prints JSON with the score and used/skipped TSS counts and exits 1 when no TSS is scored.

- The score depends on the bigWig recipe. On one real slice the same TSS set scored 11.3 with the recipe above, 14.9 with unextended read coverage and 16.0 with `--Offset 1`. Use the same recipe for every sample compared, and state it.
- `tss.bed6` needs strand in column 6 (`+` or `-`), chromosome names matching the bigWig, and the build of the BAM. A 1 bp row is used as given; for a gene interval the TSS is `start` on `+` and `end-1` on `-`. A wrong build or flipped strand makes the profile flat.
- The script is not validated against pyTSSe values. For cross-study parity use the same TSS file (for example one GENCODE protein-coding set) on both datasets.
- Pass the score to `routes/aggregate.md` as `TSS_enrichment`.

Done when `score.json` exists and the recipe and TSS file are stated.
