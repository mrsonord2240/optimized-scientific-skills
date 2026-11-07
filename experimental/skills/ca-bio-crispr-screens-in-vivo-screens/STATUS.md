# Experimental status: ca-bio-crispr-screens-in-vivo-screens

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-crispr-screens-in-vivo-screens`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Scoped to existing computational count analysis and source meta-analysis; retained GPTomics MIT attribution and computational bottleneck interpretation.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Source MLE example is saturated with one observation per animal and one shared baseline; cannot establish batch correction.
- Stouffer independence may fail with a shared control; source script clips p only below and can produce infinite z for p=1, and silently skips absent files.
- CRISPR-StAR analysis lacks a complete file-based recipe.
- Computational methods and dependencies not executed; no experimental procedures performed or added.
- Root LICENSE preflight warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-crispr-screens-in-vivo-screens). Original author:
GPTomics; upstream GPTomics/bioSkills path `crispr-screens/in-vivo-screens`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
