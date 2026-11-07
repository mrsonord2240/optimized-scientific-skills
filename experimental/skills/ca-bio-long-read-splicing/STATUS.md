# Experimental status: ca-bio-long-read-splicing

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft revises `bio-long-read-splicing`. Its earlier audit
score does not certify this changed tree. Do not present the draft as Production
Ready or Marketplace-accepted.

## Structural changes

Routed long-read workflows by intent; split alignment, microexon, discovery, curation and differential methods into topic references. Repair: made source-backed DTU and DTE routines file-based; identified unsupported custom count-only paths.

The 2026-11-07 route review replaced package-only stubs with task entry points
where the source provides a workflow. Input-dependent or incomplete source
methods are identified explicitly. New wrappers remain unexecuted.

## Validation and remaining work

The revised tree passed offline structural preflight and applicable syntax
checks. Capability coverage and command interfaces were reviewed. Tooling
preparation and meaningful behavioral validation remain before relying on results.

- Normalization only; all scientific workflows unexecuted and no behavioral audit.
- Kinnex source reports synthetic only; no real Kinnex validation.
- FLAIR ID parsing assumes gene IDs without underscores; route requires explicit confirmation.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-long-read-splicing).
Original author: GPTomics; upstream GPTomics/bioSkills path `alternative-splicing/long-read-splicing`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This modified derivative retains
the repository [MIT license](../../../LICENSE). The
[manifest](../../CANDIDATES.json) records file hashes, coverage, validation notes
and previous published identities.
