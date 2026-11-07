# Experimental status: ca-bio-crispr-screens-perturb-seq-analysis

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-crispr-screens-perturb-seq-analysis`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Routed guide assignment, Mixscape, SCEPTRE, multiome and conditional methods; retained source helpers and GPTomics MIT attribution.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Mixscape wrapper saves KO cells only, excluding controls; retain controls for downstream testing.
- Source Pertpy demo uses a built-in dataset rather than general file input; no scientific execution performed.
- Claims of unique universal FDR calibration require method and dataset review.
- Root LICENSE preflight warning: publication manifest must cite repository license evidence.
- Original source wrapper accepted a previously existing table after worker failure; repaired wrapper now refuses existing paths and requires successful worker and valid new table. Structural validation does not establish scientific correctness.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-crispr-screens-perturb-seq-analysis). Original author:
GPTomics; upstream GPTomics/bioSkills path `crispr-screens/perturb-seq-analysis`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
