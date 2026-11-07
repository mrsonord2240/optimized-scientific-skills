# Experimental status: ca-bio-proteomics-quantification

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-proteomics-quantification`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Split LFQ/MSstats/MaxLFQ, TMT bridge, SILAC and AP-MS routes; extracted source normalization and SILAC ratio workflows into file CLIs; GPTomics MIT attribution retained.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- TMT reporter extraction and impurity correction remain source snippets, not a bundled CLI.
- IRS relies on aligned protein indexes, valid nonzero reference and channel sums; verify before use.
- SILAC on/off flags describe detection, not definite biological absence.
- No scientific execution, installations or version verification.
- Root LICENSE warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-proteomics-quantification). Original author:
GPTomics; upstream GPTomics/bioSkills path `proteomics/quantification`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
