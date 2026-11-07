# Python filtering

```bash
python scripts/filter_bam.py input.bam filtered.bam -q 30 -d -p
```

This keeps mapped primary non-duplicate records at MAPQ 30; duplicates must already be marked. Use `-r chr1:1000-2000` for an indexed region and `-P` for proper pairs. For BED extraction use `python scripts/filter_by_bed.py input.bam targets.bed targets.bam`; it requires indexed coordinate-sorted BAM and BED rows of positive width. For a QNAME-consistent hash sample use `python scripts/subsample_pysam.py input.bam subset.bam 0.1 42`. See `references/pysam.md`.

Done when the output is validated and retained counts are reported.
