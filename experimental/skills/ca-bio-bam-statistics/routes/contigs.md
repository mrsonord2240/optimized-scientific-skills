# Per-contig statistics

```bash
samtools idxstats input.bam > idxstats.tsv
```

Mapped counts include secondary and supplementary records. Check contig names before mitochondrial or sex-chromosome calculations; recipes are in `references/counting-and-qc.md`.

Done when the requested contig counts or ratios are reported.
