# Duplicate rate, and filtering duplicates out

```bash
python scripts/dup_rate.py marked.bam
```

Prints total, duplicates and rate over primary alignments only (secondary and supplementary skipped). The same figure from the shell:

```bash
samtools view -c -f 1024 -F 2304 marked.bam   # duplicates
samtools view -c -F 2304 marked.bam           # primary alignments
```

`samtools flagstat` counts every record, so its percentage differs.

Without the duplicates (the marked BAM stays the master copy):

```bash
samtools view -b -F 1024 marked.bam -o nodup.bam
```

A high rate means low library complexity, over-amplification or low input. Over 50% on an assay that should not be marked at all means the wrong tool (`routes/not-markdup.md`).

Done when the rate is reported with its denominator.
