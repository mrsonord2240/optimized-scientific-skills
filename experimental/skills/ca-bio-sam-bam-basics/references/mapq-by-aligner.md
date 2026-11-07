## MAPQ Is Not Portable Across Aligners

`samtools view -q 30` does different things depending on what produced the BAM. MAPQ is an aligner-specific scale, not a universal probability:

| Aligner | MAPQ scale | "Unique" sentinel | Common gotcha |
|---------|-----------|-------------------|----------------|
| BWA-MEM / BWA-MEM2 | 0-60 | 60 | `-q 30` is sensible "high confidence" |
| minimap2 (DNA / pbmm2) | 0-60 | 60 | Spec-compliant (checked: minimap2 2.31, pbmm2 26.2.99) |
| HISAT2 | 0-60 | 60 | Spec-compliant |
| Bowtie2 | 0-42 end-to-end (0-44 with `--local`) | 42 (44 in `--local`) is the top score and most common: 97% of records in a checked run | `-q 60` drops everything; `-q 23` is a common "uniquely mapped" convention (not a probabilistic 99% threshold) |
| STAR | 0, 1, 3, 255 | **255 = uniquely mapped (sentinel, not a quality)** | `-q 255` for "unique only"; `-q 30` accidentally keeps unique only too |
| DRAGEN | 0 to `--mapq-max` (default 60) | varies | `-q 30` still meaningful; distribution shape differs (not verified here) |
| Cell Ranger / STARsolo | inherits STAR | 255 | Same trap as STAR; checked STARsolo 2.7.11b (a synthetic 10x-style run: unique reads MAPQ 255, 2-locus multimappers MAPQ 3, matching plain STAR's scale). Cell Ranger itself not run (10x Genomics gates the download behind account registration) |

MAPQ 255 means "not available" in the SAM spec; only STAR (and tools that inherit it) use it for "unique". Verify the actual scale of any unfamiliar BAM:
```bash
samtools view input.bam | awk '{print $5}' | sort -n | uniq -c | head -20   # count per MAPQ value
samtools view -H input.bam | grep '^@PG' | head -1   # which aligner produced this BAM
```
