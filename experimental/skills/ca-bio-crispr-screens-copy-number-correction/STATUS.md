# Experimental status: ca-bio-crispr-screens-copy-number-correction

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-crispr-screens-copy-number-correction`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Split diagnostic, CRISPRcleanR and panel methods; extracted diagnostic CLI, retained source example and GPTomics MIT credit.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Chronos and CERES lack complete file-based pipelines.
- CRISPRcleanR example uses KY annotation and does not export the downstream TSV named in comments.
- Source broadly calls all amplified essential hits false positives; independent validation is required.
- Methods, dependencies and API compatibility were not executed or verified.
- Root LICENSE preflight warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-crispr-screens-copy-number-correction). Original author:
GPTomics; upstream GPTomics/bioSkills path `crispr-screens/copy-number-correction`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
