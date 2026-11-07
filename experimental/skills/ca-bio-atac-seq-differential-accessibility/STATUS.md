# Experimental status: ca-bio-atac-seq-differential-accessibility

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-atac-seq-differential-accessibility`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Routed consensus analysis and conditional models, retained provider recipes and MIT notice.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Normalization only; no scientific workflows executed. Script does not apply blacklist on either branch; peak inputs must be prefiltered.
- Source spike-in snippet uses inverse spike counts as DESeq2 size factors and labels arithmetic mean as geometric; direction and scaling require validation before use.
- Source csaw and loop references are partial recipes, not a complete validated runnable analysis.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-atac-seq-differential-accessibility). Original author:
GPTomics; upstream GPTomics/bioSkills path `atac-seq/differential-accessibility`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
