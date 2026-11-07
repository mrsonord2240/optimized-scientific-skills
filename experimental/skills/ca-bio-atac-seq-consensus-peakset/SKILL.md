---
name: ca-bio-atac-seq-consensus-peakset
description: Use when building one shared peak coordinate set from per-replicate ATAC-seq peak calls before counting reads across samples.
license: MIT
category: Data Analysis
author: GPTomics
---

# Consensus Peakset

Pick the one row that matches your data. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your data and task | Read |
|--------------------|------|
| Not sure which strategy fits the question | `routes/choose.md` |
| narrowPeak per replicate, fixed-width consensus (501 bp), bedtools available | `routes/iterative-overlap.md` |
| Same, R only, no bedtools | `routes/iterative-overlap-r.md` |
| DiffBind workflow, replicate-aware fixed-width counting | `routes/diffbind.md` |
| Conditions differ strongly; keep condition-specific peaks | `routes/per-condition.md` |
| Reproducibility-filtered (IDR) set, or a quick exploratory merge | `routes/union.md` |
| A consensus exists; blacklist-filter it and make a count matrix | `routes/count.md` |
| Published peaks on another genome build | `routes/liftover.md` |
| Per-replicate peak calling, differential testing, single-cell ATAC | `routes/other-skills.md` |
| Errors or a result that looks wrong | `routes/errors.md` |

Rules on every route:

- All inputs share one genome build and chromosome naming. The blacklist and chrom.sizes match that build.
- Count every sample against the same consensus coordinates. Variable-width unions are exploratory only; they confound differential counts.
- Summit-centering needs a usable summit: narrowPeak column 10 of `-1` means none was called. The scripts refuse such rows.
- Filter against the assembly-matched blacklist before counting.
- Record the strategy, width, thresholds, genome assembly, blacklist and software versions.
- Keep the input peak files and intermediates so the set can be reproduced.

No pinned environment. featureCounts `--countReadPairs` needs Subread 2.0.2 or later.

Derived from GPTomics/bioSkills (Domen Jemec), commit d91ed3d563019e649dc854c56ccd62551359488a; MIT notice in `LICENSE`.
