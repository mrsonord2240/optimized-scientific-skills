# Duplicate Handling: pysam

## pysam Python Alternative

### Full Pipeline
`pysam.sort('-n', ...)` -> `pysam.fixmate('-m', ...)` -> `pysam.sort(...)` -> `pysam.markdup(...)` -> `pysam.index(...)`, with a record-count check, in `scripts/pysam_markdup.py`:
```bash
python scripts/pysam_markdup.py input.bam marked.bam
```

### Check Duplicate Flag
Counts primary alignments with `read.is_duplicate` (skips secondary and supplementary) and prints total, duplicates and rate, in `scripts/dup_rate.py`:
```bash
python scripts/dup_rate.py marked.bam
```

### Filter Out Duplicates
```python
import pysam

with pysam.AlignmentFile('marked.bam', 'rb') as infile:
    with pysam.AlignmentFile('nodup.bam', 'wb', header=infile.header) as outfile:
        for read in infile:
            if not read.is_duplicate:
                outfile.write(read)
```
