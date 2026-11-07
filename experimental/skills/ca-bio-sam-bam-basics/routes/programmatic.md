# Programmatic access

```bash
python scripts/view_bam.py input.bam 10
```

The arguments are file, optional limit, and optional reference FASTA; supply all three for a CRAM needing its reference. The script scans if index statistics are unavailable. For reading, writing, properties and fetch use `references/pysam.md`; region fetch requires an index. In R, Rsamtools `scanBam()` is the source-supported alternative.

Done when the requested alignments or counts are obtained with the correct coordinate convention.
