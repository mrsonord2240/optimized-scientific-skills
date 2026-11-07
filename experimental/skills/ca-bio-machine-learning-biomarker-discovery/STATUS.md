# Experimental status: ca-bio-machine-learning-biomarker-discovery

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-machine-learning-biomarker-discovery`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

File-based nested CV and stability selection using the source helper, with leakage and validation limitations.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No model training, predictions, analyses or installations.
- Stratified CV assumes independent rows; grouped and temporal data require adapted splits.
- Stability helper half-samples can fail on very rare classes and error bounds rely on assumptions; wrapper uses at least one feature budget.
- Only one null permutation is included; repeat and perform external validation.
- Skill-root LICENSE warning: cite repository license in manifest.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-machine-learning-biomarker-discovery). Original author:
GPTomics; upstream GPTomics/bioSkills path `machine-learning/biomarker-discovery`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
