# Filter flags and quality

```bash
samtools view -F 3332 -q 30 -o filtered.bam input.bam
```

Select flags and MAPQ for the assay and aligner using `references/selection.md`; the command is the original standard primary, mapped, non-duplicate MAPQ 30 filter. Verify duplicate marking first. For SV calling retain supplementary alignments and discordant pairs. FLAG choices: mapped `-F 4`, unmapped `-f 4`, proper pair `-f 2`, primary `-F 2304`, read1 `-f 64 -F 2308`, read2 `-f 128 -F 2308`, forward `-F 20`, reverse `-f 16 -F 4`. Primary is not necessarily unique. `-f` requires all bits; `-F` excludes any bit; `-G` excludes all bits set.

For complete pairs after filtering where at most two records per name remain, extract duplicated names from the filtered BAM and select with `samtools view -N`. CRAM needs `-C -T reference.fa`; SAM needs `-h`. Compare input and output counts and use `samtools quickcheck filtered.bam` and `samtools flagstat filtered.bam`.

Done when the selected output is checked and counts reported.
