# Per-contig statistics

```bash
samtools idxstats input.bam > idxstats.tsv
```

Columns are reference, length, mapped and unmapped alignment records. The mapped count includes secondary and supplementary records; for primary mapped count use `samtools view -c -F 2308 input.bam chr1`. Sum column 3 for total mapped records and sum column 4 for unmapped records. For mitochondrial fraction sum mapped records on chrM or MT and divide by the total mapped records. CRAM: use `samtools idxstats` or `pysam.idxstats`; pysam get_index_statistics can return zeros. A stale index can produce stale counts.

Done when the table and requested summaries are reported.
