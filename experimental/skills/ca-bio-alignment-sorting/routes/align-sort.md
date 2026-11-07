# Align and sort in one pass

```bash
bash scripts/align_sort.sh ref.fa R1.fq aligned.bam [R2.fq] [threads]
```

Needs `bwa index ref.fa` run once. The script pipes `bwa mem` into `samtools sort` with `pipefail`, checks that every input read appears once as a primary record (bwa can exit 0 on a truncated FASTQ), indexes, and publishes the BAM only after those checks pass. An existing output survives a failed run.

Paired: `align_sort.sh ref.fa R1.fq aligned.bam R2.fq 8`. Single-end: omit `R2.fq`.

Done when `aligned.bam` and `aligned.bam.bai` exist and flagstat printed.
