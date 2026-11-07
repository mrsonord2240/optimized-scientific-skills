# Experimental status: ca-bio-experimental-design-batch-design

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-experimental-design-batch-design`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Retained assignment and bridge scripts, added file-based SVA workflow, routed correction design context. Clarified conservative empty-cell check is not proof of non-identifiability.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific execution, installs or API verification.
- Assignment scripts reject incomplete blocks even where full model is estimable; assess rank and adapt separately.
- SVA complete-feature filtering can bias estimates; input alignment and downstream model require review.
- Skill-root LICENSE warning: cite repository license in manifest.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-experimental-design-batch-design). Original author:
GPTomics; upstream GPTomics/bioSkills path `experimental-design/batch-design`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
