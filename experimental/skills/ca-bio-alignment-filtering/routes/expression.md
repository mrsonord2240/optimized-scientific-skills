# Tags and read groups

```bash
samtools view -e '[NM] <= 5 && [AS] >= 100' -o filtered.bam input.bam
```

Select the expression or read group with `references/expressions-and-read-groups.md`. Missing tags can silently yield no records. TLEN is signed: filter both positive and negative ranges to retain both mates.

Done when the selected BAM is checked and the retained count reported.
