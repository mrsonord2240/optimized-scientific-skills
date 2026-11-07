# Experimental status: ca-bio-single-cell-cell-communication

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-single-cell-cell-communication`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Added file-based LIANA and CellChat entry points, retained CellPhoneDB and condition stability scripts with donor-level caveat.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific runs, installations or resource downloads.
- Cell-level condition permutations do not support donor-level inference; replicate-aware analysis required.
- NicheNet, spatial, Tensor-cell2cell and resource-comparison workflows need input-specific adaptation.
- LIANA source default seed and package behavior need verification.
- Skill-root LICENSE warning: cite repository license in manifest.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-single-cell-cell-communication). Original author:
GPTomics; upstream GPTomics/bioSkills path `single-cell/cell-communication`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
