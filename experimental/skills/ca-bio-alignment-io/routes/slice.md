# Slice an alignment

```bash
python scripts/alignment_io.py slice input.aln clustal output.fasta fasta --start 0 --end 100
```

Bounds are zero-based and end-exclusive. Add `--row-start` and `--row-end` to select rows. The script refuses empty or out-of-range slices.

Done when the requested region is written.
