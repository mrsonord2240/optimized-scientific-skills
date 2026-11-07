### What Each Tool Counts (and Doesn't)

| Counting category | flagstat | stats | idxstats |
|-------------------|----------|-------|----------|
| Primary alignments | `primary` line (= `in total` minus secondary minus supp) | `raw total sequences` | not separable (mapped column includes secondary/supp) |
| Secondary | `secondary` line | filtered out | counted in mapped |
| Supplementary | `supplementary` line | filtered out | counted in mapped |
| Mapping rate denominator | `in total`, *including* secondary/supp | primary only | mapped+unmapped |

flagstat lines read `<QC-passed> + <QC-failed>`. Percentages use the QC-passed column; the QC-failed reads sit in the second column and are included in `in total`, `primary` and `raw total sequences`.

For long-read data where one read produces many supplementary alignments, the senior cross-check:
```
input_read_count = flagstat_total(passed + failed) - secondary - supplementary
                 = stats_raw_total_sequences
```
Reports of "the file has 1.2M reads" where the input was actually 800k with 400k supplementary chimeric splits trace to flagstat misinterpretation.

Mate-overlap handling differs between tools, so **mean depth differs 2x on the same paired-end BAM** (test BAM: 16.77x vs 8.86x). State which convention you report:

| Tool | Overlapping mates |
|------|-------------------|
| `samtools depth`, `samtools coverage`, pysam `pileup()`, `mosdepth --fast-mode` | counted twice (16.77x); no overlap option in `coverage` |
| `samtools depth -s`, `mosdepth` (default) | counted once (8.86x) |
| `samtools mpileup` | counted once at the default `-Q 13` (8.85x); `-x` disables it (16.75x). The lower-quality mate base is zeroed, so `-Q 0` defeats the removal (16.77x). Its depth column also counts D and N positions (unlike `samtools depth`) |
| `bcftools mpileup` | counted once in `FORMAT/DP` (`-a FORMAT/DP`) and in the bases used for calling: 8.85x at the default, 16.77x with `-x` or `-Q 0`. `INFO/DP` is 16.77x either way (counted before overlap removal). Unlike `samtools mpileup`, both `FORMAT/DP` and `INFO/DP` drop to 0 across a deletion span rather than counting it (verified: a 20 bp deletion covered by 20 reads shows `samtools mpileup` depth 20 throughout, `bcftools mpileup` DP 0) |

## samtools flagstat

Fast summary of alignment flags.

```bash
samtools flagstat input.bam
```

Output:
```
10000000 + 0 in total (QC-passed reads + QC-failed reads)
9950000 + 0 primary
0 + 0 secondary
50000 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
9800000 + 0 mapped (98.00% : N/A)
9750000 + 0 primary mapped (97.99% : N/A)
9950000 + 0 paired in sequencing
4975000 + 0 read1
4975000 + 0 read2
9700000 + 0 properly paired (97.49% : N/A)
9720000 + 0 with itself and mate mapped
30000 + 0 singletons (0.30% : N/A)
15000 + 0 with mate mapped to a different chr
10000 + 0 with mate mapped to a different chr (mapQ>=5)
```
(samtools 1.13+ adds the `primary`, `primary duplicates`, and `primary mapped` lines shown above.) `properly paired` is the aligner-set proper-pair flag (both mates mapped, expected orientation and distance); `singletons` are mapped reads whose mate is unmapped; `duplicates` is non-zero only after a duplicate-marking step.

### Multi-threaded
```bash
samtools flagstat -@ 4 input.bam
```

### Machine-readable
`samtools flagstat -O tsv input.bam` prints `<passed>\t<failed>\t<label>` per line (percentages are separate `... %` rows).

### Summary table for many BAMs
```bash
printf 'Sample\tRecords\tQCfail\tPrimary\tPrimaryMapped\tProperPair\tPrimaryDup\n' > summary.tsv
for bam in *.bam; do
    sample=$(basename "$bam" .bam)
    samtools flagstat -O tsv "$bam" | awk -F'\t' -v s="$sample" '
        $3 ~ /^total/            {rec = $1 + $2; fail = $2}
        $3 == "primary"          {pri = $1 + $2}
        $3 == "primary mapped"   {pm = $1 + $2}
        $3 == "properly paired"  {pp = $1 + $2}
        $3 == "primary duplicates" {dup = $1 + $2}
        END {print s"\t"rec"\t"fail"\t"pri"\t"pm"\t"pp"\t"dup}' >> summary.tsv
done
```
`Records` is every alignment record (QC-passed + QC-failed, including secondary/supplementary); `Primary*` columns exclude secondary/supplementary and include QC-failed reads.

## samtools idxstats

Per-chromosome read counts. Uses the index; on 1.24 an unindexed BAM falls back to a slow full scan with a warning (pysam `get_index_statistics()` and mosdepth raise an error instead).

```bash
samtools idxstats input.bam
```

Output format: `chrom length mapped unmapped`
```
chr1    248956422    5000000    1000
chr2    242193529    4800000    800
chrM    16569        50000      100
*       0            0          150000
```

### Parse idxstats
Contig names differ between builds (`chrM` vs `MT`, `chrX` vs `X`), so check the names first; the recipes below exit 1 with a message instead of printing a false `0`.
```bash
# Total mapped alignments (includes secondary/supplementary)
samtools idxstats input.bam | awk '{sum += $3} END {print sum}'

# Mitochondrial percentage of mapped alignments (chrM or MT)
samtools idxstats input.bam | awk '
    $1 ~ /^(chr)?(M|MT)$/ {mt += $3; found = 1}
    {total += $3}
    END {if (!found) {print "no chrM/MT contig in idxstats" > "/dev/stderr"; exit 1}
         if (!total) {print "no mapped reads" > "/dev/stderr"; exit 1}
         printf "%.2f%% mitochondrial\n", mt/total*100}'

# Sex check (X/Y ratio; +1 avoids division by zero)
samtools idxstats input.bam | awk '
    $1 ~ /^(chr)?X$/ {x = $3; fx = 1}
    $1 ~ /^(chr)?Y$/ {y = $3; fy = 1}
    END {if (!fx || !fy) {print "no chrX/chrY contig in idxstats" > "/dev/stderr"; exit 1}
         printf "X:Y = %.2f\n", x/(y+1)}'
```

## samtools stats

Comprehensive statistics including insert size, base quality, and more.

```bash
samtools stats input.bam > stats.txt
```

### View Summary Numbers
```bash
samtools stats input.bam | grep "^SN" | cut -f 2-
```

Key summary fields:
- `raw total sequences` - Primary reads, QC-failed included (secondary/supplementary excluded)
- `reads mapped` - Mapped reads
- `reads mapped and paired` - Both mates mapped (NOT the proper-pair count)
- `reads properly paired` / `percentage of properly paired reads (%)` - Proper-pair flag
- `insert size average` / `insert size standard deviation` - Insert size mean and spread
- `inward oriented pairs` / `outward oriented pairs` / `pairs with other orientation` - Pair orientation
- `average length` - Mean read length
- `bases mapped (cigar)` - Mapped bases from the CIGAR (soft-clipped bases excluded; `bases mapped` includes them)
- `error rate` - Mismatch rate (mismatches / bases mapped (cigar))

### Generate Plots (with plot-bamstats)
```bash
samtools stats input.bam > stats.txt
plot-bamstats -p plots/ stats.txt
```

### Combine Samples (MultiQC)
```bash
multiqc . -o multiqc_out    # picks up samtools stats / flagstat / idxstats text outputs in the directory (checked on MultiQC 1.35)
```

### Stats for Specific Region
```bash
samtools stats input.bam chr1:1000000-2000000 > region_stats.txt
```

## QC Thresholds Are Assay-Specific

A single "mapping rate > 95%" rule rejects valid ATAC, ChIP, RNA-seq, metagenomics, and aDNA samples. The threshold question is "is this rate normal for this assay?" not "is this rate above 95%?" The table below is **orientation only**: approximate ranges from common practice, not sourced to a specific study and not checked on this machine's data. Judge a sample against its own assay's pipeline or facility spec.

| Metric | WGS PCR-free | WGS PCR | WES | Targeted panel | Deep panel (UMI) | RNA-seq | scRNA (10x) | ATAC | ChIP | Long-read | aDNA |
|--------|--------------|---------|-----|----------------|------------------|---------|-------------|------|------|-----------|------|
| Mapping rate | >99% | >98% | >95% | >95% | >95% | >90% | >70% | >50% | >60% | >95% | 1-50% |
| Duplicate rate | <5% | 5-15% | 20-50% | 20-50% | 50-90% pre-consensus | (skip) | (use UMI) | 10-30% | 5-30% | n/a | 20-60% |
| Proper pair rate | >95% | >95% | >85% | >80% | >80% | >70% | n/a | >50% | >70% | n/a | >60% |
| Mean MAPQ | bimodal at 0/60 | bimodal | bimodal | bimodal | bimodal | bimodal incl 255 (STAR) | 0/1/3/255 | 30-55 | 30-55 | 30-50 | 20-40 |
| Mt fraction | 0.1-2% | 0.1-2% | <1% | <0.1% | <0.1% | varies | varies | **<10% (Omni-ATAC goal; original Buenrostro-2013 libraries were often majority-mito)** | <2% | n/a | varies |

Generic red flags whatever the assay: mapping rate < 80% (contamination, wrong reference, poor quality); `error rate` > 2%.

Mean MAPQ is misleading; the distribution is bimodal (0 and aligner-max). The fraction at MAPQ >= 30 is more informative:
```bash
samtools view -c -F 2308 -q 30 in.bam   # primary, mapped, MAPQ>=30
samtools view -c -F 2308 in.bam          # primary, mapped (denominator)
# For STAR/STARsolo, use -q 255 instead of -q 30 (255 is the unique-mapping sentinel)
```
