# Subsample

```bash
samtools view -s 42.1 -b -o subset.bam input.bam
```

This keeps approximately 10% using seed 42 and selects consistently by QNAME. Use distinct seeds for sequential cuts. For target counts, `bash scripts/match_read_count.sh input.bam matched.bam --target 10000000`, or `--like normal.bam`. See `references/subsampling.md` for count definitions.

Done when the subset exists and its count is reported.
