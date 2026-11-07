# Experimental status: ca-bio-single-cell-trajectory-inference

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-single-cell-trajectory-inference`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

File-based PAGA/DPT and deterministic velocity routes plus retained CellRank script, with topology, root and model caveats.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No trajectory inference, scientific analysis, installation or API checks.
- CellRank output persistence of fate probabilities is version-dependent and must be verified on actual output.
- Palantir, Slingshot, Monocle3, tradeSeq and alternate velocity modes need study-specific adapted workflows.
- PAGA connectivity and velocity confidence are not proof of temporal direction; root and quantifier need orthogonal support.
- Skill-root LICENSE warning: cite repository license in manifest.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-single-cell-trajectory-inference). Original author:
GPTomics; upstream GPTomics/bioSkills path `single-cell/trajectory-inference`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
