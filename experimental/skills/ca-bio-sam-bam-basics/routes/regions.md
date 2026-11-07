# Retrieve regions

```bash
samtools view -M -L regions.bed -o regions.bam input.bam
```

Input must be coordinate-sorted and indexed. BED is 0-based half-open; samtools region strings are 1-based inclusive. A missing contig can produce zero records and exit 0, so check the header. Multiple overlapping region arguments duplicate records without `-M`. Python alternative: `python scripts/fetch_regions.py input.bam regions.bed > regions.sam`; add a reference as the third argument for CRAM. The Python script emits SAM records without a header.

Done when requested records are retrieved once and contig names and coordinates checked.
