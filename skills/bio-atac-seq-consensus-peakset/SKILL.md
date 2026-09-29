---
name: bio-atac-seq-consensus-peakset
description: Build a differential-ready consensus peakset from per-replicate ATAC-seq peaks using iterative overlap removal, fixed-width re-centering, and majority-rule overlap. Use when generating a stable peak coordinate system for downstream differential accessibility, ML feature engineering, cross-sample comparison, or fixed-width peak counts; covers Corces 2018 iterative overlap (501 bp), DiffBind summit re-centering, and ENCODE consistency rules.
license: MIT
---

# Consensus peakset construction

Build one shared coordinate set from per-replicate peak calls before counting reads across samples. Choose the overlap rule and peak width for the biological question; document the strategy, genome build, and filters used. Variable-width unions are exploratory and can confound differential counts.

## Workflow

1. Confirm that all input peaks use the same genome build and chromosome naming, and that summit offsets are usable when summit-centering is requested (`-1` in narrowPeak means no summit was called).
2. Choose a strategy: iterative summit-centered non-overlap for fixed-width features; DiffBind summit re-centering for its replicate-aware workflow; per-condition consensus then union when condition-specific peaks should be retained; IDR-passed union for a reproducibility-filtered set; simple merge only for exploratory work.
3. Construct the consensus, apply the selected reproducibility or overlap rule, then filter against the matching assembly blacklist before counting.
4. Count every sample against the same consensus coordinates. Convert BED to SAF only when using featureCounts; otherwise retain BED for downstream tools that accept it.
5. Record the selected method, width, thresholds, genome assembly, blacklist, and software versions.

## Routed material

- Run the Corces-style shell implementation in [`scripts/iterative_overlap.sh`](scripts/iterative_overlap.sh). It validates every selected non-comment input row as ten-column narrowPeak and rejects unusable summit offsets, including `-1`, with filename and line context before coordinate arithmetic. Review its input paths, genome sizes, and blacklist before running.
- Use [`references/method-reference.md`](references/method-reference.md) for the strategy comparison, failure modes, methods, counting examples, and citations.
- Use [`references/usage-guide.md`](references/usage-guide.md) for concise request examples.
- This Skill has no pinned environment or verified compatibility matrix. Check installed versions and command/package help before relying on version-specific options. The featureCounts `--countReadPairs` option requires Subread 2.0.2 or later.

## Inputs and outputs

Expected inputs are per-replicate peak files (narrowPeak is needed for summit-offset examples), chromosome sizes for the same assembly, and an assembly-matched blacklist when filtering. The principal output is a shared BED peakset; SAF is an optional counting annotation. Keep the original inputs and intermediate files so the final set can be reproduced.

## License and provenance

This derived Skill retains GPTomics/bioSkills material by Domen Jemec, from commit `d91ed3d563019e649dc854c56ccd62551359488a`. See [`LICENSE`](LICENSE) for the preserved MIT notice.
