# Idr



| Feature | ENCODE 3 (legacy) | ENCODE 4 (current) |
|---------|-------------------|---------------------|
| Per-rep significance threshold | `-q 0.05` directly | `-p 0.01` (loose) + IDR |
| Pseudoreplicate IDR cutoff | Not formalized | same 0.05 threshold as true replicates; rescue and self-consistency ratios |
| TSS enrichment threshold | >= 6 (older) | >= 7 (hg38, GENCODE v29) |
| Mt fraction expectation | < 25% | < 20% (Omni-ATAC < 5%) |
| Blacklist | v1 | v2 (Amemiya 2019) |
| Default genome size | hardcoded `hs`/`mm` | encouraged: deepTools effectiveGenomeSize |

To reproduce a published ENCODE 3 dataset, pin the original pipeline and threshold exactly. ENCODE 4 results are not directly numerically comparable to ENCODE 3 even on the same input BAM.

## Super-Enhancers

Super-enhancer calling (ROSE; Whyte 2013) is defined on H3K27ac, MED1 or BRD4 ChIP-seq signal, not on ATAC-seq. ROSE-style stitching of ATAC peaks can be applied when no ChIP data exist, but it is not a standard ATAC-seq analysis, and its calls are not comparable to published super-enhancer sets. For the full workflow (stitching and TSS-exclusion parameters, signal ranking and inflection, marker choice, hg38 implementation caveats), use the `bio-chipseq-super-enhancers` Skill.

## ENCODE-Style ATAC-seq Pipeline (Reference Implementation)

An ENCODE-style pattern (`scripts/call_atac_peaks.sh` runs the whole pipeline). The ENCODE pipeline instead calls peaks on Tn5-shifted tagAlign (`-f BED`) with `--call-summits`, so peak sets are comparable but not identical. `$GSIZE` is the read-length-matched effective size (see above).

```bash
# Per-replicate peak calling (loose threshold)
macs3 callpeak \
    -t rep1.filt.dedup.bam \
    -f BAM -g $GSIZE \
    -n rep1 --outdir peaks/rep1/ \
    --nomodel --shift -75 --extsize 150 \
    --keep-dup all \
    -B --SPMR \
    -p 0.01

# Pooled (all replicates)
macs3 callpeak \
    -t rep1.filt.dedup.bam rep2.filt.dedup.bam \
    -f BAM -g $GSIZE -n pooled --outdir peaks/pooled/ \
    --nomodel --shift -75 --extsize 150 --keep-dup all -B --SPMR -p 0.01

# Pseudoreplicates: disjoint halves. -s picks reads by name hash (seed.fraction), -U writes the complement,
# so no read is in both halves. Two independent -s draws with different seeds would share ~50% of reads.
samtools view -b -s 1.5 -U rep1.psr2.bam -o rep1.psr1.bam rep1.filt.dedup.bam
# (call peaks on each pseudoreplicate the same way)
```

`macs2 callpeak` takes the same flags (see `usage-guide.md` for its import failure).

`--SPMR` writes signal as Signal Per Million Reads (normalized bedGraph). `-p 0.01` is intentionally loose; IDR will tighten to a reproducible set.

## IDR for Reproducible Peaks

**Goal:** Find peaks reproducible across biological replicates at controlled IDR.

**Approach:** Rank paired peak lists by p-value, fit IDR's two-component mixture (reproducible + noise), threshold at IDR <= 0.05 for every comparison below.

```bash
# Sort peaks by p-value (column 8) so IDR scores by significance
sort -k8,8nr rep1_peaks.narrowPeak > rep1.sorted.narrowPeak
sort -k8,8nr rep2_peaks.narrowPeak > rep2.sorted.narrowPeak

# True replicates -- threshold IDR <= 0.05
idr --samples rep1.sorted.narrowPeak rep2.sorted.narrowPeak \
    --input-file-type narrowPeak --rank p.value \
    --output-file true_reps.idr \
    --idr-threshold 0.05 --plot --log-output-file idr.log

# Pseudoreplicates -- same threshold (the script runs rep1, rep2 and pooled pseudoreplicate pairs)
idr --samples psr1_peaks.narrowPeak psr2_peaks.narrowPeak \
    --input-file-type narrowPeak --rank p.value \
    --output-file psr.idr --idr-threshold 0.05 --plot
```

**ENCODE reproducibility rules** (ENCODE-DCC/atac-seq-pipeline `encode_task_reproducibility.py`). Nt = peaks passing IDR on the true replicate pair; N1, N2 = on each replicate's own pseudoreplicate pair; Np = on the pooled pseudoreplicates. Rescue ratio = max(Np, Nt) / min(Np, Nt); self-consistency ratio = max(N1, N2) / min(N1, N2). Pass if both <= 2, borderline if one > 2, fail if both > 2. IDR score >= 540 corresponds to IDR <= 0.05 (score = int(-125 log2 IDR)).

**IDR fails when:** Ranking column choice matters. `--rank p.value` (column 8) is robust; `--rank signal.value` (column 7) breaks if MACS pile-up scaling differs between replicates.
