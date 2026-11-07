# Create or refresh indices

```bash
bash scripts/ensure_index.sh input.bam
```

This checks standard index names and rebuilds when missing or any is older than the alignment, preserving existing CSI bin size. It removes existing index files during a rebuild; preserve them separately if needed. BAM defaults to BAI; for contigs longer than 2^29 bases use `samtools index -c input.bam` for CSI. CRAM uses CRAI via `samtools index input.cram`. The longest contig, not total genome length, determines the choice. A custom index name must be explicitly supplied to later commands. Sort with `samtools sort` if indexing fails on unsorted records. BAM cannot represent positions beyond 2^31-1. For multiple files pass their paths to the script. Modification times cannot detect a corrupt index or restored matching timestamps; compare idxstats to full counts when uncertain.

Done when a suitable fresh index exists and region retrieval succeeds.
