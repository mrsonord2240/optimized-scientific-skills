# Retrieve indexed regions

```bash
samtools view -o region.bam input.bam chr1:1000000-2000000
```

Regions are 1-based inclusive. Validate contig names from the header; an unknown name can return zero records with success. For multiple overlapping regions add `-M` to avoid duplicate records; for indexed BED access use `samtools view --region-file regions.bed -o region.bam input.bam`. Plain `-L` scans. For a custom index use `samtools view -X input.bam output.bai chr1:1-1000`. For CRAM supply `-T ref.fa` when the reference is not otherwise resolvable. Python alternative: `python scripts/fetch_regions.py input.bam chr1:1000000-2000000`; add `--reference ref.fa` for CRAM. This script removes stale standard indices and rebuilds them, retaining CSI type.

Done when the desired records or counts have been retrieved and checked.
