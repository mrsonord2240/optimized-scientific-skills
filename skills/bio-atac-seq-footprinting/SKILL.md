---
name: bio-atac-seq-footprinting
description: Detect transcription factor binding footprints in ATAC-seq using TOBIAS, HINT-ATAC, Wellington, or scprinter. Use when identifying bound TF sites within accessible regions, correcting Tn5 insertion bias before footprinting, choosing between cleavage-based and aggregate-based footprinters, or comparing differential TF activity between conditions.
license: MIT
category: Data Analysis
author: GPTomics
---

# TF footprinting

Detect short DNA stretches (typically 6-20 bp) of reduced Tn5 cleavage within accessible regions, where a bound TF protects DNA. This requires (1) Tn5 sequence-bias correction, (2) per-base footprint scoring, and (3) motif-anchored detection.

Tn5 has a strong sequence preference over roughly +/- 4 bp around the insertion site (Karabacak Calviello 2019). Without bias correction, "footprints" reflect enzyme preference rather than TF binding; this is the most important step.

## Workflow

1. Confirm inputs: a deduplicated, MAPQ-filtered, chrM-stripped BAM with at least about 50M nuclear reads (below that, weak or transient binders cannot be called reliably; power saturates above about 100M), a consensus peakset, a reference FASTA matching the BAM build, an assembly-matched blacklist, and a motif database (JASPAR 2024 CORE vertebrates by default; see the reference for alternatives).
2. Optionally restrict to nucleosome-free fragments (below).
3. Choose the tool for the goal: TOBIAS three-step pipeline for standard and differential two-condition work; HINT-ATAC or Wellington (`-A`) as independent second call sets; scPrinter for multi-scale scoring of a bulk library or per cell cluster. For non-model organisms use a custom motif set and a naked-DNA or retrained bias model. The decision table and tool comparison are in the reference. Per-cluster scPrinter and seq2PRINT training were only tested at low depth, where cluster footprints did not resolve (usage guide).
4. Correct Tn5 bias, score footprints, then classify bound and unbound motif sites (TOBIAS commands below).
5. Check the positive control: `scripts/run_tobias.sh` plots the CTCF aggregate per condition at all, bound, and unbound sites on uncorrected, bias-only expected (`_expected.bw`), and corrected signal. A central dip with flanking peaks at bound sites and a flatter unbound profile means the scoring separates sites sensibly; a shallow, flat, or inverted profile points to depth or bias-correction failure. Bound sites are chosen by footprint score, so their dip does not prove correction. The script therefore correlates the all-site aggregate with the bias-only expectation and exits 4 when the corrected profile still follows it (Pearson r above `BIAS_R_MAX`, default 0.2; measured values in the reference). A motif file without a CTCF motif makes the script exit 3 (`ALLOW_NO_CTCF=1` accepts that); supply another positive control then.
6. For differential work use the identical peakset and blacklist for every condition, otherwise scores are not comparable.
7. Report the tool, bias model, thresholds, genome build, motif database, and software versions. Treat single-tool calls as exploratory; for higher confidence add an independent call set (HINT-ATAC or Wellington footprints, or ChIP-seq peaks) and report site-level concordance with `scripts/site_concordance.sh`, always beside the unbound-site control (definition and measured values in the reference). Absence of a footprint is not absence of binding, especially for nuclear receptors, pioneer factors, and short-motif TFs; see the per-TF caveats in the reference.

## Tn5 cut geometry

Tn5 cuts with a 9 bp stagger. Apply a +4 bp shift on the + strand and a -5 bp shift on the - strand to read 5' ends before per-base counting. TOBIAS and HINT-ATAC do this internally, scPrinter applies it according to its `--shift` setting, and custom counting must apply it (deepTools `alignmentSieve --ATACshift` applies the same shift). Skipping it produces roughly 9 bp asymmetry in aggregate footprints.

## TOBIAS three-step pipeline

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

## scPrinter multi-scale footprints (GPU advised)

`scripts/scprinter_footprint.py` runs the tested classic route: genome Tn5 bias prediction, fragment import, and `get_footprint_score` at scales 2-100 around each region centre, for one pseudobulk or, with `--groups barcode<TAB>group`, per cell cluster. `--shift` is the shift already applied to the fragment ends: `0,0` (default) for raw BAM-derived fragments, `4,-5` for Cell Ranger fragments. Each `mode` is a footprint scale in bp (footprint and flank radius default to it); in the tested bulk chr1 slice CTCF separated bound from unbound sites at modes 10-30 but not at 50. Missing inputs exit 2. Setup pins, the fragment-file recipe, per-cluster limits and seq2PRINT training are in the usage guide. Compare bound versus unbound motif sites, as with TOBIAS.

## Nucleosome-free fragment filtering

Sub-100 bp fragments carry most TF signal. Filtering sharpens footprints but discards nucleosome-borne information, so keep the unfiltered BAM for other analyses.

```bash
samtools view -h sample.bam | \
    awk 'substr($0,1,1)=="@" || ($9 > 0 && $9 < 100) || ($9 < 0 && $9 > -100)' | \
    samtools view -b > sample.nfr.bam
samtools index sample.nfr.bam
```

## Routed material

- Run the two-condition TOBIAS pipeline with per-condition CTCF positive-control plots, the bias-correction check (exit 4), and a ranked differential summary (|change| among p <= 0.05, `PMAX` to change) through [`scripts/run_tobias.sh`](scripts/run_tobias.sh). Its positional arguments are the two BAMs, peaks, genome FASTA, blacklist, motif file, output directory, and cores. Review input paths and the motif file name before running.
- Compare two call sets with [`scripts/site_concordance.sh`](scripts/site_concordance.sh) (bound sites, unbound sites, footprints BED).
- Run bulk or per-cluster scPrinter scoring with [`scripts/scprinter_footprint.py`](scripts/scprinter_footprint.py).
- Use [`references/method-reference.md`](references/method-reference.md) for the tested-with versions, tool taxonomy, Tn5 bias and alternatives, per-TF failure modes, goal decision table, call-set concordance, motif databases, in silico variant effects, common errors, citations, and related Skills.
- Use [`references/usage-guide.md`](references/usage-guide.md) for per-tool environments, RGT and scPrinter setup, fragment preparation, per-cluster scPrinter, seq2PRINT training, and request examples.

## License and provenance

This derived Skill retains GPTomics/bioSkills material by Domen Jemec, from commit `d91ed3d563019e649dc854c56ccd62551359488a`. See [`LICENSE`](LICENSE) for the preserved MIT notice.
