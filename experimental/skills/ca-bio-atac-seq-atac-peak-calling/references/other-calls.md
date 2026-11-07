# Other-Calls



| Scenario | Recommended caller | Why |
|----------|-------------------|-----|
| Bulk ATAC, 2-3 reps, depth >= 25M | MACS3/MACS2 ENCODE-style pipeline + IDR | Reproducible, comparable to published peaksets |
| Bulk ATAC, 1 sample (no rep) | MACS3 callpeak with `-q 0.05`; do not run IDR | IDR is meaningless without reps; tighter q-value substitutes |
| Bulk ATAC, deep library, want NFR + flanking nuc structure | MACS3 hmmratac | HMM separates NFR from nucleosome flanks |
| Multi-replicate joint analysis where rep weighting is symmetric | Genrich `-j` ATAC mode | Joint p-value across reps; built-in chrM and blacklist |
| Cell type with broad accessible domains | MACS3 `--broad --broad-cutoff 0.1` for domains; narrow for differential | Domain-level inference vs site-level; super-enhancer calling is a separate ChIP-based workflow (see Super-Enhancers) |
| FFPE / degraded chromatin (flat fragment dist) | MACS3 callpeak with stringent `-q 0.01`; never HMMRATAC | HMM needs fragment periodicity |
| scATAC pseudobulk per cluster | MACS3 callpeak per cluster + iterative overlap | See bio-atac-seq-single-cell-atac |
| Want fixed-width consensus peaks for differential | Call broadly, then re-center to summits +/- 250 bp | See bio-atac-seq-consensus-peakset |
| Plant / non-model organism | MACS3 with `-g <effective_size>`; verify size empirically | Default `-g hs/mm` invalid; compute via khmer |

## Reconciliation: When Callers Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| MACS narrow peaks much fewer than Genrich | Different statistical models and peak merging (GM12878 chr1 slice: Genrich 14,328 joint peaks at `-q 0.05` vs MACS3 4,803 for one replicate; `-q 0.01` gave 22,936, so lowering `-q` did not reduce the count) | Do not tune `-q` for parity; compare by overlap with IDR peaks or keep the MACS-supported subset |
| HMMRATAC misses peaks MACS finds | Library too shallow OR fragment-size periodicity weak | Trust MACS; HMMRATAC is depth-sensitive |
| HMMRATAC calls peaks MACS misses | HMM is sensitive to mid-strength accessibility flanked by phased nucleosomes | Inspect; often genuine but unconfirmed by short-fragment signal |
| Same peak called by all but width 2x different | Broad mode vs narrow mode mismatch | Standardize: re-center to summit +/- 250 bp for differential |
| Per-rep MACS calls peak; pooled MACS does not | One rep dominates; lambda smoothes it out in pooled | Trust pooled + IDR over per-rep counts |

**Operational rule for high-confidence reporting:** Require a peak to pass IDR <= 0.05 on true replicates AND survive blacklist/greylist filtering AND have mean signalValue >= 5 across reps. Two callers from different families (MACS + Genrich) agreeing within 250 bp is acceptable evidence when IDR is unavailable.

## Blacklist and Greylist

```bash
# ENCODE blacklist (Amemiya 2019) -- always remove
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
bedtools intersect -v -a peaks.narrowPeak -b hg38-blacklist.v2.bed > peaks.no_blacklist.narrowPeak

# Sample-specific greylist (input-derived high-signal regions; rarely available for ATAC)
# For ATAC, ENCODE recommends pooling all samples' top-percentile signal and removing
# regions exceeding 100x median coverage as a "soft greylist"
```

Blacklist is mandatory; greylist is optional and most useful when the same library prep produces consistent artifact regions across samples.

## NFR-Only Peak Calling

**Goal:** Call peaks using only sub-nucleosomal fragments (<100 bp) for sharper TF-binding-relevant accessibility.

**Approach:** Pre-filter BAM to short fragments, then call peaks with parameters scaled to the smaller fragment length.

```bash
samtools view -h sample.dedup.bam | \
    awk 'substr($0,1,1)=="@" || ($9 > 0 && $9 < 100) || ($9 < 0 && $9 > -100)' | \
    samtools view -b > nfr.bam
samtools index nfr.bam

macs3 callpeak -t nfr.bam -f BAM -g $GSIZE -n sample_nfr \
    --nomodel --shift -37 --extsize 75 \
    --keep-dup all -p 0.01
```

`--shift -37 --extsize 75` halves both parameters to match shorter fragments; this is a fragment-scaled convention for NFR-focused input, not a TOBIAS-specified setting (TOBIAS instead applies the +4/-5 Tn5 correction to the full BAM via ATACorrect).

## Output Files (narrowPeak)

| Column | Field | Notes |
|--------|-------|-------|
| 1-3 | chrom, start, end | 0-based, half-open |
| 4 | name | MACS auto-numbers |
| 5 | score | Min(int(-10*log10(qvalue)), 1000) |
| 6 | strand | `.` for ATAC |
| 7 | signalValue | Fold enrichment over local lambda |
| 8 | pValue | -log10 p |
| 9 | qValue | -log10 q (BH-FDR) |
| 10 | summit_offset | Peak summit relative to start |

Convert to bigWig for browsers: `sort -k1,1 -k2,2n sample_treat_pileup.bdg > sample.sorted.bdg && bedGraphToBigWig sample.sorted.bdg chrom.sizes sample.bw` (bedGraphToBigWig is multi-pass and cannot read from a pipe/stdin, so sort to a file first).

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `--shift/--extsize ignored` warning | Used `-f BAMPE` with these flags | Switch to `-f BAM` or remove the flags |
| 0 peaks called | Forgot `--nomodel`; MACS tries to build a shifting model and fails | Add `--nomodel --shift -75 --extsize 150` |
| Peak count >> 500k | Did not deduplicate; or did not remove chrM; or `-q` too loose | Pre-filter (samtools view -F 1804 -q 30; samtools idxstats); use `-q 0.01` |
| `Sequence chrM not found` (Genrich) | Wrong chromosome name in `-e` flag (chrM vs MT) | Match BAM header naming convention |
| HMMRATAC out of memory | Deep library / large genome; the current tool is `macs3 hmmratac` (Python, not Java) | Increase available RAM and use a scratch `--outdir`; the `-Xmx`/`HMMRATAC.jar` heap flags apply only to the deprecated standalone Java HMMRATAC |
| Peaks shifted by 75 bp from expected positions | Forgot `--shift -75` (cuts at one end of read) | Add the shift; positions are now centered on Tn5 cut site |
| IDR returns 0 reproducible peaks | Sorted by wrong column; ranks are random | Sort each peakset by `-k8,8nr` (p-value descending) |
