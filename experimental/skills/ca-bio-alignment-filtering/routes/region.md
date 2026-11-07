# Region or BED

```bash
samtools view -L targets.bed -o targets.bam input.bam
```

BED coordinates are 0-based half-open. A region argument such as `chr1:1000-2000` is 1-based inclusive and requires an index; check contigs using `samtools idxstats` because an absent name can yield zero reads with success. Add `-P` and an index to fetch mates outside the selected region. Combine with `-q` or FLAG options as appropriate.

Done when the output is validated and its count compared with the input.
