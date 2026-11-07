# Experimental status: ca-bio-single-cell-splicing

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft revises `bio-single-cell-splicing`. Its earlier audit
score does not certify this changed tree. Do not present the draft as Production
Ready or Marketplace-accepted.

## Structural changes

Replaced package-only route stubs with file-based task entry points; source-incomplete paths remain explicitly identified.

The 2026-11-07 route review replaced package-only stubs with task entry points
where the source provides a workflow. Input-dependent or incomplete source
methods are identified explicitly. New wrappers remain unexecuted.

## Validation and remaining work

The revised tree passed offline structural preflight and applicable syntax
checks. Capability coverage and command interfaces were reviewed. Tooling
preparation and meaningful behavioral validation remain before relying on results.

- All three newly added wrappers are syntax-only and require tooling validation; no scientific workflows executed.
- MARVEL modality and donor-aware group comparisons require adapting the source reference.
- Sierra, BRIE quantification and pseudobulk aggregation retain multi-step reference workflows; SpliZ source workflow was previously unrun.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-single-cell-splicing).
Original author: GPTomics; upstream GPTomics/bioSkills path `alternative-splicing/single-cell-splicing`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This modified derivative retains
the repository [MIT license](../../../LICENSE). The
[manifest](../../CANDIDATES.json) records file hashes, coverage, validation notes
and previous published identities.
