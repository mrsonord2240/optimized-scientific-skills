# FLAIR differential splicing

```bash
bedtools bamtobed -bed12 -i aligned.bam > aligned.bed
```

Continue correct, collapse, quantify and diffSplice using `references/flair.md`. Provide orthogonal junctions for novel splice sites. For each tested event type verify its DRIMSeq output exists; exit zero alone is insufficient.

Done when requested event counts and tests are available and reviewed.
