# Inputs



`ScoreBigwig` was formerly `FootprintScores`. `ATACorrect` writes `<bam-prefix>_uncorrected.bw`, `_bias.bw`, `_expected.bw`, and `_corrected.bw` per input BAM; substitute the real prefix below.

```bash
# Step 1: bias correction, per condition, same peaks and blacklist
TOBIAS ATACorrect \
    --bam cond1.bam --genome hg38.fa \
    --peaks consensus.bed --blacklist hg38-blacklist.v2.bed \
    --outdir cond1_corrected/ --cores 16

# Step 2: continuous per-base footprint score
TOBIAS ScoreBigwig \
    --signal cond1_corrected/cond1_corrected.bw \
    --regions consensus.bed \
    --output cond1_footprints.bw \
    --cores 16

# Step 3: motif-anchored bound/unbound calls and differential
TOBIAS BINDetect \
    --motifs JASPAR2024_CORE_vertebrates.pfm \
    --signals cond1_footprints.bw cond2_footprints.bw \
    --genome hg38.fa --peaks consensus.bed \
    --outdir bindetect/ \
    --cond-names cond1 cond2 \
    --cores 16
```

BINDetect reports per-motif bound counts, mean scores, `cond1_cond2_change`, and `cond1_cond2_pvalue`. A positive change with a low p-value means more bound in cond1; a negative change means more bound in cond2. The change is a difference in mean footprint score across motif sites, not a fold-change; there is no formal cutoff, so calibrate against positive controls (observed magnitudes are in the reference). Further interpretation patterns are in the reference.

## Independent second call sets

Both tools take the ATAC BAM directly (paired-end) and need a regions BED; HINT-ATAC also needs RGT data set up once (usage guide).

```bash
rgt-hint footprinting --atac-seq --paired-end --organism=hg38 \
    --output-location hint_out --output-prefix sample reads.bam peaks.bed   # hint_out/sample.bed

wellington_footprints.py -A peaks.bed reads.bam wellington_out/             # -A = ATAC mode; *FDR*.bed files
```



Sub-100 bp fragments carry most TF signal. Filtering sharpens footprints but discards nucleosome-borne information, so keep the unfiltered BAM for other analyses.

```bash
samtools view -h sample.bam | \
    awk 'substr($0,1,1)=="@" || ($9 > 0 && $9 < 100) || ($9 < 0 && $9 > -100)' | \
    samtools view -b > sample.nfr.bam
samtools index sample.nfr.bam
```
