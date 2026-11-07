# DNA BAM in Python (pysam)

```bash
python scripts/pysam_markdup.py in.bam marked.bam
```

Sorts by name, runs `fixmate -m`, sorts by coordinate, runs `markdup`, indexes, and exits with an error if the record count changed. It has no assay gate: check `routes/not-markdup.md` first if the assay is not plain DNA. pysam has no Windows wheel; run on Linux, macOS or WSL.

Do not re-sort or rewrite the BAM in Python between `fixmate` and `markdup`: the `ms` and `MC` tags are lost and samtools stops with an error.

Filter duplicates out of a marked BAM:

```python
import pysam
with pysam.AlignmentFile('marked.bam', 'rb') as infile, \
     pysam.AlignmentFile('nodup.bam', 'wb', header=infile.header) as outfile:
    for read in infile:
        if not read.is_duplicate:
            outfile.write(read)
```

Done when `marked.bam` and its index exist and the script printed matching record counts.
