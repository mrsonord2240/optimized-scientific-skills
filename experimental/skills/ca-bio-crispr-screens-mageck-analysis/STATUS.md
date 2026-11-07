# Experimental status: ca-bio-crispr-screens-mageck-analysis

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-crispr-screens-mageck-analysis`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Split existing read count, RRA, MLE and visualization workflows; extracted source time-course and volcano functions with file interfaces and retained GPTomics MIT attribution.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Time-course source function computes same-sign flags but does not require them; route warns.
- Source example shell uses hardcoded FDR positions and is not copied; outputs need header-based interpretation.
- MAGeCKFlute and VISPR are source snippets not bundled runnable wrappers; package compatibility remains unverified.
- MLE permutation FDR is stochastic and design rank/confounding require review.
- No scientific execution or installation.
- Root LICENSE preflight warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-crispr-screens-mageck-analysis). Original author:
GPTomics; upstream GPTomics/bioSkills path `crispr-screens/mageck-analysis`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
