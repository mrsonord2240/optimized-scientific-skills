# Read counts and QC

```bash
samtools flagstat input.bam > flagstat.txt
```

For insert size, base quality and error summaries run `samtools stats input.bam > stats.txt`. Flagstat totals include secondary and supplementary records; primary excludes them. See `references/counting-and-qc.md` for multi-sample tables, assay-specific expectations and reporting conventions, and `references/qc-pitfalls.md` for blind spots.

Done when the requested metrics and their denominators are reported.
