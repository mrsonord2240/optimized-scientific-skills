---
name: bio-atac-seq-atac-peak-calling
description: Call accessible chromatin regions from ATAC-seq BAM files using MACS3, MACS2, Genrich, or HMMRATAC. Use when identifying open chromatin from aligned ATAC-seq, choosing between point-source vs HMM peak callers, applying ENCODE-style pseudoreplicate IDR, removing blacklist regions, or re-centering peaks on summits for downstream differential analysis.
license: MIT
tool_type: cli
primary_tool: macs3
category: Data Analysis
author: GPTomics
---

# ATAC-seq Peak Calling

Call accessible chromatin regions from ATAC-seq BAM files. Identify Tn5-hypersensitive open chromatin, treating fragments as point insertion events (not protein-bound regions as in ChIP-seq) and accounting for the lack of input control.

## Workflow

1. Verify upstream BAM is deduplicated, MAPQ-filtered, and chrM-removed.
2. Choose a caller based on depth, replicate count, and downstream use: MACS2/MACS3 for ENCODE-style point-source calls; Genrich for joint-replicate mode; MACS3 hmmratac for HMM-based fragment-class separation; HOMER for convenience in downstream motif analysis.
3. Decide on shift-extend vs BAMPE fragment modeling: `-f BAM --shift -75 --extsize 150` (ENCODE pattern, point-source) vs `-f BAMPE` (paired-fragment, biologically-scaled widths).
4. Set effective genome size from deepTools table, not the legacy `-g hs/mm` shorthand.
5. For the ENCODE-style pipeline: call per-replicate and pooled peaks, split each replicate into disjoint pseudoreplicates, run IDR on true replicates and pseudoreplicates at one threshold (0.05), and judge reproducibility by the rescue and self-consistency ratios (rule in `method-reference.md`).
6. Filter against ENCODE blacklist (Amemiya 2019 v2) and optionally a sample-derived greylist.
7. Produce narrowPeak + optional bigWig signal track outputs. The script writes the conservative (true-replicate IDR) set only.

## Routed material

- Run the ENCODE-style reference implementation in [`scripts/call_atac_peaks.sh`](scripts/call_atac_peaks.sh). It checks that the BAMs are readable, indexed and free of chrM reads, then runs per-replicate and pooled MACS calls, disjoint pseudoreplicates, IDR, the rescue and self-consistency ratios, and blacklist filtering of the conservative set. Genome size and blacklist are required arguments; `MACS=macs2` selects MACS2. Its usage header lists the ENCODE differences.
- Use [`references/method-reference.md`](references/method-reference.md) for algorithmic taxonomy, failure modes, caller decision tree, reconciliation patterns, ENCODE vs ENCODE 3 differences, and citations.
- Use [`references/usage-guide.md`](references/usage-guide.md) for quick-start prompts, caller quick reference, and effective genome size tables.
- Install with the two-environment recipe in `usage-guide.md` (IDR needs numpy<1.24, MACS3 needs newer). `--version` and `--help` can succeed on a broken MACS2, so confirm a real `callpeak` run. Tested on GM12878 chr1 data with MACS3 3.0.4, MACS2 2.2.9.1 (import fails on glibc >= 2.31 without an `*_finite` LD_PRELOAD shim), Genrich 0.6.2, samtools 1.24, bedtools 2.31.1, IDR 2.0.4.2, bedGraphToBigWig 482. ROSE, HOMER, chromap and standalone Java HMMRATAC were not run.

## Inputs and outputs

Expected inputs are BAM files (deduplicated, chrM-removed, MAPQ-filtered), chromosome sizes for the assembly, and an assembly-matched blacklist (ENCODE v2 mandatory, greylist optional). The principal output is narrowPeak peaks; bigWig signal tracks and IDR scores are optional. Keep original inputs and intermediate files so final peak sets can be reproduced.

## License and provenance

This derived Skill retains GPTomics/bioSkills material by Domen Jemec, from commit `d91ed3d563019e649dc854c56ccd62551359488a`. See [`LICENSE`](LICENSE) for the preserved MIT notice.
