# pysam

```python
import pysam
pysam.sort('-o', 'sorted.bam', 'input.bam')
pysam.sort('-N', '-o', 'namesorted.bam', 'input.bam')   # -n for samtools-only consumers
pysam.sort('-@', '4', '-m', '2G', '-T', '/tmp/sortpfx', '-o', 'sorted.bam', 'input.bam')
pysam.merge('-c', '-p', '-f', 'merged.bam', 'a.bam', 'b.bam', 'c.bam')
```

- Do not load records into a list and call `sorted()`; it exhausts memory around 30M reads (about 10 GB human BAM). `pysam.sort` uses samtools' external merge sort.
- Check order with `scripts/check_order.py`; re-sort only when it fails (`routes/check-order.md`).
- Merge only inputs sorted the same way, with RG IDs unique (`routes/merge.md`).
- Piping from an aligner needs `pipefail`, or a crashed aligner still returns 0 and leaves a truncated BAM:

```python
import subprocess
subprocess.run(['bash', '-o', 'pipefail', '-c',
                'bwa mem ref.fa reads.fq | samtools sort -o aligned.bam'], check=True)
```

Done when `samtools quickcheck -v` passes on the output.
