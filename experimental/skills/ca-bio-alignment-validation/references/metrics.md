## Proper Pairing Rate

Percentage of mapped paired reads that are correctly paired. Every count below is over primary records (`-F 2304`; mapped: `-F 2308`), the same denominators as `flagstat` "primary".

### samtools flagstat

```bash
samtools flagstat input.bam

samtools flagstat input.bam | grep "properly paired"
```

### Calculate Pairing Rate

```bash
proper=$(samtools view -c -f 2 -F 2308 input.bam)
paired=$(samtools view -c -f 1 -F 2308 input.bam)
awk -v p="$proper" -v n="$paired" 'BEGIN{if (n) printf "Proper pairing rate: %.2f%%\n", 100*p/n; else print "no paired reads"}'
```

Multiply before dividing: `bc` with `scale=2` on `$proper / $paired * 100` truncates the quotient first (99.96% prints as 99.00%).

## Strand Balance

The metric is the **forward fraction F/(F+R)** over mapped primary reads. A balanced 0.48-0.52 applies to WGS / WES / generic DNA-seq on autosomes (do not use the F/R quotient, which is ~1.0 when balanced). Expected to deviate for: stranded RNA-seq (deliberately strand-asymmetric -- verify with RSeQC `infer_experiment.py -r genes.bed12 -i rna.bam`, needs a BED12 gene model: ~0.5/0.5 = unstranded, one fraction near 1 = stranded; RSeQC 5.0.4), bisulfite (CT vs GA), small-RNA / strand-specific RNA-seq, and chrY/chrM regions. Per-chromosome strand imbalance >5% on autosomes is a field-convention rule of thumb (no single primary citation) — it picks up aligner artifacts; on chrX/chrY it suggests sex-mismatch.

### Calculate Strand Fraction

```bash
forward=$(samtools view -c -F 2324 input.bam)          # mapped, primary, forward
reverse=$(samtools view -c -f 16 -F 2308 input.bam)    # mapped, primary, reverse
echo "Forward: $forward"
echo "Reverse: $reverse"
awk -v f="$forward" -v r="$reverse" 'BEGIN{if (f+r) printf "Forward fraction: %.3f\n", f/(f+r)}'
```

### Check Strand Bias per Chromosome

Needs an index (`samtools index input.bam`).


```bash
samtools idxstats input.bam | awk '$1!="*" && $3>0 {print $1}' | head -25 | while read -r chr; do
    fwd=$(samtools view -c -F 2324 input.bam "$chr")
    rev=$(samtools view -c -f 16 -F 2308 input.bam "$chr")
    awk -v c="$chr" -v f="$fwd" -v r="$rev" 'BEGIN{if (f+r) printf "%s: F=%d R=%d forward fraction=%.3f\n", c, f, r, f/(f+r)}'
done
```

## Mapping Quality Distribution

Mapped primary reads only (`-F 2308`): unmapped records carry MAPQ 0 (or a stray value) and would distort both the histogram and the mean.

### Extract MAPQ Distribution

```bash
samtools view -F 2308 input.bam | cut -f5 | sort -n | uniq -c | sort -k2 -n
```

### Calculate Mean MAPQ

```bash
samtools view -F 2308 input.bam | awk '{sum+=$5; count++} END {if (count) print "Mean MAPQ:", sum/count}'
```

### MAPQ Distribution Is Bimodal and Aligner-Specific

Mean MAPQ is misleading; distributions are bimodal (0 and aligner-max). For aligner-specific scales and "unique mapping" sentinels, see bio-sam-bam-basics. The fraction of primary mapped reads at MAPQ >= 30 (`samtools view -c -F 2308 -q 30` over `-c -F 2308`) is a more informative summary than the mean.

## Chromosome Coverage Balance

### Calculate Per-Chromosome Coverage

Mapped reads per bp (a density, not depth):
```bash
samtools idxstats input.bam | awk '$2>0 {printf "%s\t%.4f\n", $1, $3/$2}' | head -25
```

### Check for Aneuploidy / Sex Chromosome Imbalance

Median-normalized per-autosome density (1.0 = expected diploid; 0.5 = monosomy/sex; 1.5 = trisomy). Plain POSIX awk (no gawk `asort`), plain `chr1`-`chr22` / `1`-`22` only (alt, decoy, HLA and unplaced contigs would distort the median):
```bash
samtools idxstats in.bam | awk '$2>0 && $3>0 && $1 ~ /^(chr)?[0-9]+$/ {print $1, $3/$2}' \
  | sort -k2,2g \
  | awk '{c[NR]=$1; v[NR]=$2} END {
      if (NR==0) { print "no autosomes with reads" > "/dev/stderr"; exit 1 }
      med = v[int(NR/2)+1]
      for (i=1; i<=NR; i++) printf "%s\t%.3f\n", c[i], v[i]/med }'
```

For full ancestry / contamination / relatedness checking, use `verifybamid2` or `somalier` -- they account for population AFs, not just per-contig depth.

## Mismatch Rate

### Picard CollectAlignmentSummaryMetrics

```bash
picard CollectAlignmentSummaryMetrics \
    I=input.bam \
    R=reference.fa \
    O=alignment_summary.txt
```

### Key Metrics

Picard reports these as fractions (0.9996 = 99.96%); the Good values are percentages, so multiply by 100.

| Metric | Description | Good Value |
|--------|-------------|------------|
| PCT_PF_READS_ALIGNED | Mapped fraction | > 95% (0.95) |
| PF_MISMATCH_RATE | Mismatches | < 1% (0.01) |
| PF_INDEL_RATE | Indels | < 0.1% (0.001) |
| STRAND_BALANCE | Strand ratio | ~0.5 |

## Quality Thresholds Summary

Germline short-read DNA defaults; the example validators grade mapping rate, proper pairing, strand balance and mean MAPQ with these bands. Per-assay values (RNA-seq, ATAC, ChIP, long-read, aDNA, ...) are in bio-bam-statistics.

| Metric | Good | Warning | Fail |
|--------|------|---------|------|
| Mapping rate | > 95% | 90-95% | < 90% |
| Proper pairing | > 90% | 80-90% | < 80% |
| Singletons (`flagstat`) | < 5% | 5-10% | > 10% |
| Duplicate rate (assay-specific) | see bio-bam-statistics decision table | -- | -- |
| Strand balance (forward fraction) | 0.48-0.52 | 0.45-0.55 | Outside |
| Mean MAPQ | > 40 | 30-40 | < 30 |
| GC bias (`NORMALIZED_COVERAGE`) | < 1.2x | 1.2-1.5x | > 1.5x (field-convention bands; Picard CollectGcBiasMetrics does not prescribe specific cutoffs) |

## When a Metric Fails

| Symptom | Check |
|---------|-------|
| Low mapping rate (< 90%) | Correct reference genome version; sample species matches reference; contamination with other species; read quality/adapters before alignment |
| Poor proper pairing (< 80%) | Insert size vs aligner expectation; chimeric reads from structural variants; library protocol; re-align with adjusted parameters |
| Abnormal insert size | Bimodal = mixed libraries; very wide = degraded DNA; adapter contamination; compare with the library table |
| Strand imbalance | Capture bias in targeted sequencing; PCR artifacts; look at specific chromosomes for local bias |
| GC bias | Strong correlation = PCR amplification issues, library complexity or low-input DNA; may need GC correction before analysis |
| Low MAPQ | Repetitive / duplicated regions with multi-mappers; longer reads or paired-end help; filter low MAPQ for variant calling |
