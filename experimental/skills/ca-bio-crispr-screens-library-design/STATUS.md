# Experimental status: ca-bio-crispr-screens-library-design

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-crispr-screens-library-design`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Routed existing TSS window and generic oligo formatting helpers; retained design and QC references and GPTomics MIT attribution. Synthetic demo is explicitly identified as incomplete.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Source candidate example uses synthetic sequences and controls; no real input-file pipeline or off-target scorer is supplied.
- GC heuristic is not Rule Set 2; source reports Azimuth compatibility problems.
- Oligo helper does not validate sequence or change output based on vector option.
- No new organism-specific procedures or scientific executions performed.
- Root LICENSE preflight warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-crispr-screens-library-design). Original author:
GPTomics; upstream GPTomics/bioSkills path `crispr-screens/library-design`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
