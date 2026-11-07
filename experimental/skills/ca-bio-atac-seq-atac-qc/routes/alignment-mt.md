# Alignment rate and mitochondrial fraction

```bash
samtools flagstat sample.bam
samtools idxstats sample.bam | awk '{m+=$3} $1=="chrM"{c=$3} END{print "chrM fraction of mapped:", c/m}'
```

Shown, not scripted. The mt fraction is chrM mapped reads over all mapped reads. If the BAM has no `chrM`, the mitochondrial contig is named `MT`, `Mt` or similar: match the name from `idxstats` to your filter.

- Above 50% chrM: re-prep (Omni-ATAC pre-spin or digitonin lysis with mt removal). If 30-50%, the library may be salvageable by removing chrM but yield is reduced.
- Report the value; feed it to `routes/aggregate.md` as `mt_fraction`.

Done when alignment rate and mt fraction are recorded.
