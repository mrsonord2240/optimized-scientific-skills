# Duplicate Handling: UMI-aware deduplication

## UMI-Aware Deduplication

For UMI libraries (10x scRNA, ctDNA panels, Twist/IDT/Roche UMI capture), naive markdup destroys information. Use UMI-aware tools:

### umi_tools dedup

Input must be **coordinate-sorted and indexed**.
Pass `--paired` for paired-end libraries: without it the mates are deduplicated independently and the output is silently wrong (5689 vs 2805 records on a paired-end capture BAM). Add `--random-seed=1` for a reproducible output: umi_tools picks among tied reads at random, so without it the count moves by about 1 between runs (5688 or 5689 here).

`scripts/umi_tools_dedup.sh` runs both forms with `--method=directional --random-seed=1` and fails on an empty output:
```bash
# 10x / scRNA: group by cell barcode (CB) + UMI (UB). With absent CB/UB, --per-cell writes an EMPTY BAM and
# still exits 0, so the script checks for CB:Z: in the first 1000 records and exits 2 if there are none.
bash scripts/umi_tools_dedup.sh scrna cellranger_possorted.bam dedup.bam

# Bulk UMI, paired-end (UMI in the RX tag): sorts and indexes a copy first, passes --paired
bash scripts/umi_tools_dedup.sh bulk-paired raw.bam dedup.bam
```

### fgbio consensus (bulk UMI / ctDNA, best practice for low-VAF detection)

`GroupReadsByUmi` needs the mate mapping-quality (`MQ`) tag on every read (see Common Errors). `samtools fixmate -m` on name-grouped input adds it; alternatively `fgbio SetMateInformation` on queryname-sorted input. Consensus reads are written **unmapped**; re-align them before variant calling. Single-strand and duplex use different grouping strategies and are separate branches. `scripts/fgbio_consensus.sh` name-sorts, adds the mate tags, groups and calls consensus:

```bash
# Single-strand molecular consensus (--strategy=adjacency --edits=1, then CallMolecularConsensusReads --min-reads=1)
bash scripts/fgbio_consensus.sh single raw.bam consensus.bam

# Duplex (xGen-Prism, NEBNext duplex): --strategy=paired, then CallDuplexConsensusReads --min-reads 1 1 0
bash scripts/fgbio_consensus.sh duplex raw.bam duplex.bam
```

Duplex needs `--strategy=paired`, which requires RX as two UMIs joined by `-` (UMI1-UMI2; a single-UMI RX fails with IllegalArgumentException, so single-UMI libraries use the single-strand branch) and writes MI tags with `/A` `/B` strand suffixes. `CallDuplexConsensusReads` on adjacency-grouped reads crashes (StringIndexOutOfBoundsException). If the UMI is in a separate FASTQ instead of the RX tag, annotate first and pass the annotated BAM: `fgbio AnnotateBamWithUmis -i raw.bam -f umi.fastq -o annotated.bam`.

### Picard UMI-aware marking
```bash
picard UmiAwareMarkDuplicatesWithMateCigar I=coordsort_fixmate.bam O=marked.bam M=metrics.txt \
    UMI_METRICS=umi_metrics.txt UMI_TAG_NAME=RX
```

`--method=directional` is the default and correct -- do not use `--method=unique`, which treats single-base UMI errors as different molecules. `samtools markdup --barcode-tag RX` (UMI/barcode handling added in samtools 1.16) does exact-match UMI grouping; adequate for IDT xGen Duplex but insufficient for single-UMI applications where 1-edit errors are common.
