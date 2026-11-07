# UMI libraries

Plain `samtools markdup` destroys UMI data. Pick the branch, then run its command.

| Library | Command |
|---------|---------|
| 10x or scRNA, cell barcode CB and UMI UB tags | `bash scripts/umi_tools_dedup.sh scrna in.bam out.bam` |
| Bulk paired-end, UMI in the RX tag, deduplicate | `bash scripts/umi_tools_dedup.sh bulk-paired in.bam out.bam` |
| Single-strand molecular consensus (ctDNA, UMI capture) | `bash scripts/fgbio_consensus.sh single in.bam out.bam` |
| Duplex consensus (xGen-Prism, NEBNext duplex) | `bash scripts/fgbio_consensus.sh duplex in.bam out.bam` |

- `umi_tools_dedup.sh` uses `--method=directional --random-seed=1` and fails on an empty output. `scrna` exits 2 if there is no `CB:Z:` tag in the first 1000 records (otherwise `--per-cell` writes an empty BAM and exits 0). `bulk-paired` sorts and indexes a copy and passes `--paired`; without it mates are deduplicated separately and the output is silently wrong. Never `--method=unique`.
- `fgbio_consensus.sh` name-sorts, runs `samtools fixmate -m` (the MQ tag fgbio needs), groups and calls consensus. Consensus reads are written unmapped: re-align them before variant calling.
- Duplex needs RX as two UMIs joined by `-` (UMI1-UMI2). A single-UMI RX must use `single`. The mode picks the grouping strategy because the wrong one crashes.
- UMI in a separate FASTQ: `fgbio AnnotateBamWithUmis -i raw.bam -f umi.fastq -o annotated.bam`, then run the script on the annotated BAM.
- Picard alternative for Twist/IDT/Roche UMI capture: `picard UmiAwareMarkDuplicatesWithMateCigar I=coordsort_fixmate.bam O=marked.bam M=metrics.txt UMI_METRICS=umi_metrics.txt UMI_TAG_NAME=RX`. It is BETA and silently keeps duplicates on transcriptome-aligned BAMs; do not use it for RNA.
- `samtools markdup --barcode-tag RX` matches UMIs exactly: acceptable for xGen Duplex, not for single-UMI libraries where one-base UMI errors are common.

Done when the output BAM is non-empty and the script printed its record count.
