# View and count

```bash
samtools view -h input.bam | head -100
```

Use `samtools view -H input.bam` for header only. `samtools view -c input.bam` counts alignment records, including secondary and supplementary; add `-F 2304` for primary alignments. `-h` retains the header when writing SAM. Use `-@ N` for additional threads. For CRAM reference requirements see `references/cram-reference.md`.

Done when the requested records, header or counts are reported.
