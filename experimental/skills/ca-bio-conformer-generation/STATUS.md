# Experimental status: ca-bio-conformer-generation

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft revises `bio-conformer-generation`. Its earlier audit
score does not certify this changed tree. Do not present the draft as Production
Ready or Marketplace-accepted.

## Structural changes

Separated conformer, macrocycle, CREST and property workflows; preserved source implementation and default settings. Repair: replaced import-only task routes with file-based workflows or a documented unsupported boundary.

The 2026-11-07 route review replaced package-only stubs with task entry points
where the source provides a workflow. Input-dependent or incomplete source
methods are identified explicitly. New wrappers remain unexecuted.

## Validation and remaining work

The revised tree passed offline structural preflight and applicable syntax
checks. Capability coverage and command interfaces were reviewed. Tooling
preparation and meaningful behavioral validation remain before relying on results.

- Normalization only; no scientific execution or behavioral audit.
- Source CREST platform and CLI claims are version-specific; confirm installed documentation before use.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-conformer-generation).
Original author: GPTomics; upstream GPTomics/bioSkills path `chemoinformatics/conformer-generation`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This modified derivative retains
the repository [MIT license](../../../LICENSE). The
[manifest](../../CANDIDATES.json) records file hashes, coverage, validation notes
and previous published identities.
