# Fragment-size periodicity, TSSEscore, FRiP

```bash
picard CollectInsertSizeMetrics I=filtered.bam O=isize.txt H=isize.pdf M=0.5
Rscript scripts/atac_qc_metrics.R sample.bam peaks.narrowPeak sample_qc TxDb.Hsapiens.UCSC.hg38.knownGene
```

Arguments: BAM, peaks (narrowPeak or BED; `-` skips FRiP), output prefix, TxDb package (default hg38 knownGene). It needs a coordinate-sorted, duplicate-marked, paired-end BAM, and a TxDb whose seqlevel style matches the BAM (it stops otherwise; use the TxDb of the BAM's genome). It counts only MAPQ >= 30, non-duplicate, properly paired, non-chrM alignments, and writes `<prefix>_fragsize.pdf` and `<prefix>_metrics.csv` (nuclear reads, TSSEscore, NFR and mono fractions, periodicity class, FRiP, grades).

- Periodicity classes are `3+ peaks`, `NFR+mono`, `flat/single`: a smoothed-density heuristic, not ENCODE-defined. Read the PDF for finer patterns (`routes/diagnose.md`).
- ATACseqQC `TSSEscore` (100 bp centre, 1000 bp flanks) is not comparable to ENCODE-style TSS enrichment and is typically 2-3x larger. For the ENCODE score use `routes/tss-encode.md`.
- Insert size is valid only for paired-end data; for single-end switch to deepTools `fragmentSize`.
- FRiP is the fraction of read pairs overlapping a peak; report the peak caller, version and cutoff.

Done when `<prefix>_metrics.csv` exists.
