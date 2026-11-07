## Picard CollectRnaSeqMetrics and Gene-Body Coverage

**Goal:** Mapping distribution (coding / UTR / intronic / intergenic / rRNA) and 5'-3' bias.

Picard needs a `refFlat` file. Build it from the BED12 (checked with Picard 3.5.0):

```bash
awk 'BEGIN{OFS="\t"} {n=$10; split($11,sz,","); split($12,st,","); s=""; e="";
     for(i=1;i<=n;i++){s=s ($2+st[i]) ","; e=e ($2+st[i]+sz[i]) ","}
     print $4,$4,$1,$6,$2,$3,$7,$8,n,s,e}' genes.bed12 > refFlat.txt

picard CollectRnaSeqMetrics \
    I=sample.bam O=sample.rna_metrics.txt REF_FLAT=refFlat.txt \
    STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
    RIBOSOMAL_INTERVALS=rRNA_intervals.interval_list      # header-bearing interval_list, not BED

geneBody_coverage.py -i sample.bam -r genes.bed12 -o sample_geneBody --skip-plot
```

Strand flags per library (foot-gun):

| Library | rMATS `--libType` | featureCounts | Picard `STRAND_SPECIFICITY` |
|---------|-------------------|---------------|------------------------------|
| Reverse-stranded (TruSeq Stranded, NEB Ultra II Directional; dUTP) | `fr-firststrand` | `-s 2` | `SECOND_READ_TRANSCRIPTION_STRAND` |
| Forward-stranded (Lexogen QuantSeq FWD, some ligation kits) | `fr-secondstrand` | `-s 1` | `FIRST_READ_TRANSCRIPTION_STRAND` |
| Unstranded | `fr-unstranded` | `-s 0` | `NONE` |

On a planted dUTP library `SECOND_READ_TRANSCRIPTION_STRAND` gave `PCT_CORRECT_STRAND_READS` 1.0 and `FIRST_READ...` 0.0. STAR has no library-strand flag; use `--outSAMstrandField intronMotif` for `XS` tags.

| Metric | Healthy (convention) | Concerning |
|--------|----------------------|------------|
| PCT_CODING_BASES | >=50% | <30% (degradation or mis-priming) |
| PCT_UTR_BASES | 20-40% | >50% (3' bias) |
| PCT_INTRONIC_BASES | <30% (poly(A)); <60% (rRNA-depleted) | >50% in poly(A): pre-mRNA contamination |
| PCT_INTERGENIC_BASES | <10% | >20%: genomic DNA contamination |
| PCT_RIBOSOMAL_BASES | see rRNA section | |
| MEDIAN_5PRIME_TO_3PRIME_BIAS | 0.7-1.3 | **<0.5 = possible 3' bias** (degraded RNA or poly(A) capture); **>2 = possible 5' bias** |
| Gene-body coverage curve | Flat | Strong 3' skew = low RIN or library mis-prep |

Checked on planted BAMs: a 3'-biased library gave `MEDIAN_5PRIME_TO_3PRIME_BIAS` 0.047 (uniform 1.009) and `geneBody_coverage` 5'/3' 0.04 vs 1.02. With 3' bias junction reads concentrate at the 3' end and miss internal junctions. Treat the Picard ratio as a screen, not a standalone call: on a low-depth chrX subset it was 0.2745 while the gene-body 5'/3' ratio was 0.847. Confirm a flagged ratio with gene-body coverage and inspect the number of usable reads/expressed transcripts before declaring a bias.

## Strandedness Verification

```bash
infer_experiment.py -i sample.bam -r genes.bed12 -s 200000     # BED6 or BED12
```

Read the two `Fraction of reads explained by ...` lines (and `failed to determine`):

| Output | Library | rMATS `--libType` |
|--------|---------|-------------------|
| ~0.5 / ~0.5 | Unstranded | `fr-unstranded` |
| >=0.9 `"1++,1--,2+-,2-+"` (PE) or `"++,--"` (SE) | Forward-stranded | `fr-secondstrand` |
| >=0.9 `"1+-,1-+,2++,2--"` (PE) or `"+-,-+"` (SE) | Reverse-stranded (dUTP / TruSeq stranded) | `fr-firststrand` |

A wrong strand setting halves usable junction reads, so verify before quantification. Leaky libraries: 20% and 40% planted leakage read 0.80 and 0.59 on the dominant string. At 0.7-0.9 report the leakage and expect that fraction of reads on the wrong strand; below 0.7 treat the library as unstranded.

## Annotation Choice

| GENCODE level | Contents | Use for |
|---------------|----------|---------|
| Basic | High-confidence canonical isoforms | Standard rMATS, leafcutter, SUPPA2; event-level AS |
| Comprehensive | All transcripts including putative/predicted | Transcript-level DTU (DRIMSeq+DEXSeq, satuRn), isoform discovery |
| RefSeq | NCBI curated | Less complete than GENCODE |
| Ensembl | Same content as GENCODE in vertebrates | Different attribute conventions |

Comprehensive captures rare isoforms but adds annotation noise and multiple-testing burden; basic can under-detect rare isoforms in DTU.

## rRNA Contamination Check

**Post-alignment (fraction of primary mapped records overlapping rRNA):**

```bash
total=$(samtools view -c -F 0x904 sample.bam)                       # -F 0x904: drop unmapped, secondary, supplementary
rrna=$(samtools view -c -F 0x904 -L rRNA_intervals.bed sample.bam)
awk -v r=$rrna -v t=$total 'BEGIN{printf "rRNA: %.1f%%\n", 100*r/t}'
```

Without `-F 0x904` the denominator includes secondary and unmapped records: a planted 22.2% library read 15.6% and passed the 20% rule. Picard `PCT_RIBOSOMAL_BASES` (0.222 on the same BAM) is the alternative.

**Pre-alignment:** `fastq_screen --aligner minimap2 --conf fastq_screen.conf --threads 8 sample_R1.fq.gz` (bowtie2, bowtie and bwa work too if installed). Setup facts, verified with fastq_screen 0.16.0:
- The default aligner is bowtie2; without it (or another aligner) the command exits 255.
- `DATABASE<TAB>rRNA<TAB>/path/to/prefix` gives the index basename. For minimap2 fastq_screen opens `prefix.fa.gz` (gzipped FASTA next to the `.mmi`); if it is missing the run exits 0 with only an `Aligner warning` in the log and reports 100% unmapped.
- Read `sample_R1_screen.txt`: the `%One_hit_one_genome` plus `%Multiple_hits_*` columns of the rRNA row are the rRNA fraction (25.00% on a planted 25% library).

| rRNA fraction | Library | Status |
|---------------|---------|--------|
| <5% | rRNA-depleted | Excellent |
| 5-20% | rRNA-depleted | Acceptable; some leakage |
| >=20% | rRNA-depleted | Failed depletion; redo |
| <5% | poly(A) | Healthy |
| >5% | poly(A) | Suggests degraded RNA or poor selection |
