# FASTA index

```bash
samtools faidx reference.fa
```

Retrieve a 1-based inclusive interval with `samtools faidx reference.fa chr1:1000-2000`. In pysam, `FastaFile.fetch('chr1',999,2000)` uses 0-based half-open coordinates.

Done when the FAI exists and any requested interval is retrieved.
