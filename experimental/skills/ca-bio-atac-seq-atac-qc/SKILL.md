---
name: ca-bio-atac-seq-atac-qc
description: Use when judging whether an ATAC-seq library passes ENCODE 4 quality criteria or deciding which replicates to drop before peak calling.
tool_type: mixed
primary_tool: deeptools
license: MIT
category: Data Analysis
author: GPTomics
---

# ATAC-seq Quality Control

Pick the one row that matches your data. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your data and task | Read |
|--------------------|------|
| Which metrics to compute and what passes (ENCODE 4 thresholds) | `routes/thresholds.md` |
| Alignment rate and mitochondrial fraction from a BAM | `routes/alignment-mt.md` |
| NRF, PBC1, PBC2 from a duplicate-retained BAM | `routes/complexity.md` |
| Fragment-size periodicity, ATACseqQC TSSEscore, FRiP from a paired-end BAM | `routes/fragment-tss-frip.md` |
| ENCODE-style TSS enrichment from a bigWig | `routes/tss-encode.md` |
| Replicate correlation and fingerprints | `routes/replicates.md` |
| Would re-sequencing rescue a low-complexity library (preseq) | `routes/preseq.md` |
| Grade the metrics, one report card per sample, MultiQC | `routes/aggregate.md` |
| A metric failed and you need the cause; flat fragment sizes; odd scores | `routes/diagnose.md` |
| Sex check, spike-in normalization, cell-cycle effect | `routes/extra-checks.md` |
| Peak calling, nucleosome positioning, single-cell QC, upstream FASTQ QC | `routes/other-skills.md` |

Rules on every route:

- Compute NRF and PBC on the raw mapped BAM before deduplication (the script applies MAPQ >= 30 and drops chrM). A deduplicated BAM gives NRF 1.0 trivially.
- State which TSS enrichment implementation and bigWig recipe produced a score. ENCODE-style and ATACseqQC `TSSEscore` differ by 2-3x, and the bigWig recipe shifts the score too. Use one recipe for every sample compared.
- Match the mitochondrial contig name (`chrM`, `MT`, `M`) to `samtools idxstats` before filtering.
- FRiP depends on the peak set: pin the caller version and cutoff before comparing.
- The ENCODE mt-fraction bands assume Omni-ATAC; standard ATAC often shows 30-50% chrM. Above 50% is the practical reject line.
- NRF under 0.9 alone is not fatal (legitimate duplicates at hyperaccessible sites); confirm bottlenecking with PBC1 < 0.7, PBC2 < 1.0 and coverage pile-ups.
- Do not re-run QC on a chrM-stripped BAM to hide a high mt fraction.
- ENCODE thresholds are organism-specific; for non-model organisms use cohort percentile rank.
- Correlate replicates with Spearman, not Pearson, with `--skipZeros`.

Tested with deepTools 3.5, Picard 3.1, samtools 1.19, bedtools 2.31, ATACseqQC 1.26, pysam 0.22, pyBigWig 0.3, numpy 1.26, pandas 2.2, MultiQC 1.21.

Derived from GPTomics/bioSkills, commit d91ed3d563019e649dc854c56ccd62551359488a; MIT notice in `LICENSE`.
