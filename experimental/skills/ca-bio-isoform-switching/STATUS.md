# Experimental status: ca-bio-isoform-switching

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft revises `bio-isoform-switching`. Its earlier audit
score does not certify this changed tree. Do not present the draft as Production
Ready or Marketplace-accepted.

## Structural changes

Split workflow by statistical intent and functional annotation; retained original executable scripts and their defaults. Repair: made source-backed DTU and DTE routines file-based; identified unsupported custom count-only paths.

The 2026-11-07 route review replaced package-only stubs with task entry points
where the source provides a workflow. Input-dependent or incomplete source
methods are identified explicitly. New wrappers remain unexecuted.

## Validation and remaining work

The revised tree passed offline structural preflight and applicable syntax
checks. Capability coverage and command interfaces were reviewed. Tooling
preparation and meaningful behavioral validation remain before relying on results.

- Normalization only; scripts and scientific commands unexecuted; no behavioral audit.
- Raw count route and small-sample null behavior remain source-described limitations.
- Count-only IsoformSwitchAnalyzeR and single-cell import are source-incomplete or untested; route marks custom or sibling workflow.

- Licensed annotators and single-cell import remain unverified.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-isoform-switching).
Original author: GPTomics; upstream GPTomics/bioSkills path `alternative-splicing/isoform-switching`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This modified derivative retains
the repository [MIT license](../../../LICENSE). The
[manifest](../../CANDIDATES.json) records file hashes, coverage, validation notes
and previous published identities.
