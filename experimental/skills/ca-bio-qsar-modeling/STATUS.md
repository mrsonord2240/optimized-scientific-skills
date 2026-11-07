# Experimental status: ca-bio-qsar-modeling

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-qsar-modeling`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Intent routes for training, validation, uncertainty, calibration and interpretation; retained source command defaults.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Normalization only; no scientific execution or behavioral audit.
- Chemprop and MAPIE commands are source-version-specific; confirm installed interfaces before use.
- GP active learning and graph attribution require adaptation to user data; no packaged graph attribution.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-qsar-modeling). Original author:
GPTomics; upstream GPTomics/bioSkills path `chemoinformatics/qsar-modeling`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
