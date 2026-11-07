# Experimental status: ca-bio-crispr-screens-base-editing-analysis

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-crispr-screens-base-editing-analysis`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Routed existing post-screen analysis, predictions and annotation; retained source scripts and moved conditional references.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Normalization only; no scientific execution or behavioral audit.
- Spacer helper can return rows without a target position; route requires filtering.
- Efficiency helper skips missing sample tables and classifies any loss of target nucleotide as editing; inspect assay schema.
- BE-Hive and CRISPResso interfaces are source-date-specific.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-crispr-screens-base-editing-analysis). Original author:
GPTomics; upstream GPTomics/bioSkills path `crispr-screens/base-editing-analysis`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
