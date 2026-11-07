# Experimental status: ca-bio-machine-learning-atlas-mapping

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-machine-learning-atlas-mapping`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Parameterized source scANVI reference/query template and retained marker check, including empty-verdict guard.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No model training, prediction, download, installation or scientific analysis.
- Neighbor disagreement and group marker checks cannot guarantee out-of-distribution detection or label identity.
- Reference overlap, held-out validation and model training need study-specific review.
- Other mappers and saved-reference-only surgery require separate adaptation.
- Skill-root LICENSE warning: cite repository license in manifest.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-machine-learning-atlas-mapping). Original author:
GPTomics; upstream GPTomics/bioSkills path `machine-learning/atlas-mapping`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
