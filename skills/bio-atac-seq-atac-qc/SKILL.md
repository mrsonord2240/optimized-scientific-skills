---
name: bio-atac-seq-atac-qc
description: ATAC-seq library quality control -- TSS enrichment, FRiP, fragment-size periodicity, library complexity (NRF/PBC1/PBC2), mitochondrial fraction, and ENCODE 4 thresholds. Use when assessing whether an ATAC-seq library passes ENCODE acceptance criteria, diagnosing transposition artefacts, comparing Omni-ATAC vs standard prep quality, or selecting which replicates to drop before peak calling.
tool_type: mixed
primary_tool: deeptools
license: MIT
category: Data Analysis
author: GPTomics
---

## Version Compatibility

Reference examples tested with: deepTools 3.5+, Picard 3.1+, samtools 1.19+, bedtools 2.31+, ATACseqQC 1.26+, pysam 0.22+, pyBigWig 0.3+, numpy 1.26+, pandas 2.2+, MultiQC 1.21+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt.

# ATAC-seq Quality Control

**"Does my ATAC library pass ENCODE quality criteria?"** -> Compute the seven canonical metrics (depth, alignment rate, mitochondrial fraction, library complexity, fragment-size periodicity, TSS enrichment, FRiP) and compare against ENCODE 4 thresholds, then diagnose failures.

- CLI: `picard CollectInsertSizeMetrics`, `samtools flagstat`, `samtools idxstats`
- CLI: `deeptools plotFingerprint`, `computeMatrix reference-point` + `plotProfile`
- R: `ATACseqQC::TSSEscore`, `ATACseqQC::fragSizeDist`, `ATACseqQC::PTscore`
- Python: fragment-level NRF/PBC (pysam); pyBigWig for TSS enrichment

## Workflow

1. Collect alignment statistics: `samtools flagstat` for total alignment rate, `samtools idxstats` for chrM fraction.
2. Calculate library complexity: run [`scripts/library_complexity.py`](scripts/library_complexity.py) on the duplicate-retained BAM. It counts distinct fragments (paired-end; single-end reads by 5' position) after MAPQ >= 30 and chrM exclusion; see the script docstring for options.
3. Analyze fragment sizes: `picard CollectInsertSizeMetrics I=filtered.bam O=isize.txt H=isize.pdf M=0.5`, then run [`scripts/atac_qc_metrics.R`](scripts/atac_qc_metrics.R) `<bam> [peaks.narrowPeak|-] [prefix] [TxDb package]` for the periodicity class, ATACseqQC TSSEscore and FRiP. It needs a duplicate-marked paired-end BAM and a TxDb matching the BAM's genome (default hg38).
4. Calculate TSS enrichment: build the bigWig (below), then run [`scripts/encode_tss_enrichment.py`](scripts/encode_tss_enrichment.py) `<bw> <tss.bed6>`. Score depends on the bigWig recipe; it exits 1 when no TSS is scored.
5. Aggregate metrics: run [`scripts/aggregate_qc.py`](scripts/aggregate_qc.py) on a metrics JSON to write a sample-wide `*_mqc.tsv` with per-metric and overall grades (missing metrics grade `NA`, overall `INCOMPLETE`).
6. Compare replicates: use deepTools `multiBamSummary` + `plotCorrelation` for Spearman correlation; `plotFingerprint` for enrichment visualization.

See [`references/usage-guide.md`](references/usage-guide.md) for quick-start prompts and examples.

## ENCODE 4 ATAC-seq Acceptance Thresholds

| Metric | Definition | Ideal | Acceptable | Reject | Source |
|--------|-----------|-------|------------|--------|--------|
| Nuclear reads (after dedup, no chrM) | Mapped reads (2 per fragment), MAPQ >= 30, non-chrM, deduped | >= 50M | 25-50M | < 25M | ENCODE 4 minimum: 25M single-end reads, 50M paired-end reads (25M fragments); tiers are this Skill's convention |
| Alignment rate | Mapped / total reads | >= 95% | 80-95% | < 80% | ENCODE 4 |
| Mitochondrial fraction | chrM / total mapped | < 5% (Omni-ATAC), < 20% (standard) | 20-50% | > 50% | Working convention (Corces 2017); ENCODE 4 sets no mt threshold |
| NRF (Non-Redundant Fraction) | Distinct fragments / total fragments (single-end: 5' positions / reads) | >= 0.9 | 0.7-0.9 | < 0.7 | ENCODE 4 ideal > 0.9; other bands are working convention (Landt 2012 defines the metric) |
| PBC1 (PCR Bottlenecking Coefficient 1) | Fragments seen once / distinct fragments | >= 0.9 | 0.7-0.9 | < 0.7 | ENCODE 4 ideal > 0.9; other bands are working convention |
| PBC2 | Fragments seen once / fragments seen twice | >= 3.0 | 1.0-3.0 | < 1.0 | ENCODE 4 ideal > 3; other bands are working convention |
| TSS enrichment (hg38) | Avg signal at TSS / avg flanking | >= 7 | 5-7 | < 5 | ENCODE 4 (cutoffs vary by annotation; ENCODE lists GRCh38 RefSeq) |
| FRiP (Fraction Reads in Peaks) | Reads in MACS peaks / total | >= 0.3 | 0.2-0.3 | < 0.2 | ENCODE 4 |
| Insert-size periodicity | NFR + mono-nuc + di-nuc peaks visible | Clear 3+ peaks | NFR + mono only | Flat / single peak | ENCODE 4 requires NFR and mononucleosome peaks; Buenrostro 2013 |

ENCODE thresholds are organism-specific. Mouse (mm10, GENCODE M21) TSS enrichment >= 5 is acceptable; non-model organisms have no published threshold (use cohort percentile rank instead). Methodology evolves; verify against the current ENCODE ATAC-seq Standards before reporting.

## TSS Enrichment: ENCODE Method vs ATACseqQC Method

The two most common implementations DO NOT produce identical scores.

| Method | Numerator | Denominator | Scaling |
|--------|-----------|-------------|---------|
| ENCODE-style (`encode_tss_enrichment.py`; not validated against pyTSSe-reported values) | Mean signal in 100 bp window centered at TSS | Mean signal in 100 bp window at +/- 1900 to +/- 2000 bp (flanks) | Per-base normalization to flanks; reported as fold-enrichment |
| ATACseqQC TSSEscore | Sum signal in TSS +/- 100 bp | Sum signal at +/- 1000 bp flanking windows | Different window sizes; ratios are larger |
| deeptools plotProfile | Visual; numeric ratio not standardized | Reference-point matrix | No standard score; for visualization only |

**Trigger:** Comparing a TSS score across studies.

**Mechanism:** Different normalization windows shift the absolute number; ATACseqQC's TSSEscore is typically 2-3x ENCODE's because of the wider flank.

**Symptom:** Reported score 21 vs ENCODE-ideal 7 mismatch. Likely the calculator was ATACseqQC; the equivalent ENCODE score might be 8.

**Fix:** State which implementation and bigWig recipe were used. For ENCODE-style scoring use [`scripts/encode_tss_enrichment.py`](scripts/encode_tss_enrichment.py) (or Kundaje-lab pyTSSe).

**bigWig recipe:** the score is recipe-dependent. On one real slice the same TSS set scored 11.3 with the recipe below, 14.9 with unextended read coverage and 16.0 with `--Offset 1`. Use one recipe for every sample being compared:

```bash
bamCoverage -b filtered.bam -o sample.bw -bs 1 --extendReads --normalizeUsing None
```

The TSS BED must be BED6 with strand in column 6 (1 bp TSS rows are used as given; for gene intervals the TSS is `start` on `+`, `end-1` on `-`), with chromosome names matching the bigWig.

## Fragment-Size Periodicity Patterns

`scripts/atac_qc_metrics.R` reports only three classes (`3+ peaks`, `NFR+mono`, `flat/single`; smoothed-density heuristic, not ENCODE-defined). The finer patterns below are for reading the plot.

| Pattern | Visual signature | Interpretation | Action |
|---------|-----------------|----------------|--------|
| Strong tri-modal | NFR (~50bp) >> mono (~200bp) > di (~400bp) > tri (~600bp) peaks | Excellent transposition; well-positioned chromatin | Pass |
| Clear bi-modal | NFR + mono only, di and tri faint | Acceptable; common in Omni-ATAC | Pass |
| Single broad peak | Flat after NFR or no NFR | Over-transposition (too much Tn5) OR degraded chromatin | Reject; cannot distinguish nucleosomes |
| Inverted (mono >> NFR) | Mono peak dominant, NFR weak | Under-transposition OR chromatin condensation | Caution; peak counts will be low |
| Sharp 147 bp spike with no flanks | Tight peak at 147 bp | ChIP-seq input contamination (MNase-like) | Reject; not ATAC-grade |
| 10.4 bp helical periodicity overlay | Sub-peaks at 50, 60, 70, 80 bp on NFR | Excellent chromatin structure resolution; helical phasing visible | Pass; high-quality |

The 10.4 bp helical periodicity is a Buenrostro 2013 hallmark: it reflects the helical pitch of B-form DNA, with Tn5 preferring outward-facing minor grooves on nucleosomal DNA. Its presence is a positive QC indicator but not required.

## Per-Metric Failure Modes and Diagnostics

See [`references/method-reference.md`](references/method-reference.md) for detailed failure mode analysis, mitochondrial fraction issues, library complexity bottlenecking, TSS enrichment troubleshooting, FRiP interpretation, replicate correlation failures, sex-chromosome QC, cell-cycle effects, spike-in normalization, and common error resolution.

## Cross-Replicate QC

```bash
# Spearman correlation (more robust than Pearson for ATAC)
multiBamSummary bins -bs 10000 -p 8 \
    --bamfiles rep1.bam rep2.bam rep3.bam \
    -o multi.npz

plotCorrelation -in multi.npz \
    --corMethod spearman --whatToPlot heatmap --skipZeros \
    -o spearman_heatmap.png

# Fingerprint (per-bin signal cumulative -- diagonal = no enrichment, sharp curve = good)
plotFingerprint -p 8 -b rep1.bam rep2.bam rep3.bam \
    --labels rep1 rep2 rep3 \
    --skipZeros --numberOfSamples 50000 \
    -o fingerprint.png \
    --outQualityMetrics fingerprint_metrics.txt
```

deepTools fingerprint quality metrics report a synthetic JS distance without a reference; the (non-synthetic) Jensen-Shannon distance column is only computed when a reference sample is supplied via `--JSDsample`. Larger values indicate stronger enrichment.

## Library Complexity Extrapolation (preseq)

**Goal:** Predict whether re-sequencing would rescue a low-NRF library, separating "library is bottlenecked" from "we just sequenced too shallow."

**Approach:** Fit preseq's rational-function (Pade) approximation of the Good-Toulmin power-series estimator on observed BAM positions; extrapolate distinct-fragment yield as a function of additional sequencing depth. Input must be a coordinate-sorted, duplicate-retained BAM (a deduplicated BAM has no redundancy to extrapolate). Paired-end ATAC needs `-P` (fragments); without it preseq counts mates as independent reads.

```bash
# c_curve: observed complexity at current depth; step -s must be well below the fragment count
# (a step above depth returns only the "0 0" row), e.g. ~1/10 of it
preseq c_curve -B -P -s 1e5 -o sample.ccurve.tsv sample.bam

# lc_extrap: predicted complexity at higher depth (-e here 200M; preseq default -e is 1e10, step -s default 1M)
preseq lc_extrap -B -P -e 200000000 -s 5000000 -o sample.lcextrap.tsv sample.bam
```

Interpretation: if `lc_extrap` shows distinct-fragment count flattening before 100M reads, the library is bottlenecked (re-sequencing won't help; re-prep needed). If it continues to climb, re-sequencing will recover more unique reads. Use alongside NRF/PBC1/PBC2 to decide library re-prep vs deeper sequencing.

## MultiQC Aggregation

```bash
# Run after generating per-sample QC outputs
multiqc \
    fastqc/ \
    picard/ \
    samtools_stats/ \
    macs2/ \
    deeptools/ \
    -o multiqc_report
```

MultiQC ingests Picard CollectInsertSizeMetrics, samtools flagstat, deepTools plotFingerprint output, and MACS peaks tables. It does NOT compute TSS enrichment or NRF; pipe a custom `_mqc.tsv` for those via [`scripts/aggregate_qc.py`](scripts/aggregate_qc.py).

## Common Errors

See [`references/method-reference.md`](references/method-reference.md) for error diagnosis table, TSS enrichment scoring mismatches, NRF/PBC computation errors, mitochondrial fraction naming issues, insert-size flatness, and FRiP drift solutions.

## References

- Buenrostro JD et al 2013 Nat Methods 10:1213 (ATAC-seq protocol; fragment-size periodicity)
- Corces MR et al 2017 Nat Methods 14:959 (Omni-ATAC; mt fraction reduction protocol)
- Landt SG et al 2012 Genome Res 22:1813 (ENCODE/modENCODE QC framework, NRF/PBC definitions; the PBC1/PBC2 split is a later ENCODE-pipeline refinement)
- ENCODE 4 ATAC-seq Data Standards (encodeproject.org/atac-seq) -- canonical thresholds
- Ou J et al 2018 BMC Genomics 19:169 (ATACseqQC R package; TSSEscore implementation)
- Ramirez F et al 2016 Nucleic Acids Res 44:W160 (deepTools, plotFingerprint JSD)
- Daley T & Smith AD 2013 Nat Methods 10:325 (preseq library-complexity extrapolation model; the lc_extrap re-sequencing decision)

## Related Skills

- atac-seq/atac-peak-calling - FRiP requires peaks; QC drives accept/reject before calling
- atac-seq/nucleosome-positioning - Fragment-size analysis
- atac-seq/single-cell-atac - per-cell QC has different thresholds
- read-qc/quality-reports - upstream FastQC
- alignment-files/bam-statistics - samtools flagstat / idxstats
- alignment-files/duplicate-handling - dedup before NRF/PBC computation

## License and provenance

This derived Skill retains GPTomics/bioSkills material from commit `d91ed3d563019e649dc854c56ccd62551359488a`. See [`LICENSE`](LICENSE) for the preserved MIT notice.
