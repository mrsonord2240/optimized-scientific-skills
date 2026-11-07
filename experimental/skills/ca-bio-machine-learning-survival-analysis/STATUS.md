# Experimental status: ca-bio-machine-learning-survival-analysis

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-machine-learning-survival-analysis`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Added file-based fixed-penalty Cox research fit and Aalen-Johansen CIF, preserving prediction validation and censoring caveats.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No training, prediction, analyses, installations or API checks.
- Cox command fits fixed penalty only; it does not tune, validate, calibrate or produce an individual clinical predictor.
- Fine-Gray, RSF, nested survival validation and calibration need study-specific workflows.
- Aalen-Johansen tie handling may jitter and must be reviewed.
- Skill-root LICENSE warning: cite repository license in manifest.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-machine-learning-survival-analysis). Original author:
GPTomics; upstream GPTomics/bioSkills path `machine-learning/survival-analysis`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
