# Experimental status: ca-bio-crispr-screens-combinatorial-screens

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-crispr-screens-combinatorial-screens`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Split scoring, interaction-model and architecture guidance into routes; retained GPTomics MIT credit and source method caveats.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Library design, counting, Perturb-seq and annotation lack complete source implementations; custom work or canonical sibling skills required.
- Source makes conflicting claims about MLE interaction FDR and formal testing; route flags uncertainty.
- Z-score null and BH interpretation require scientific review; source uses fixed filenames and overwrites duplicate singleton IDs.
- Scientific execution, dependencies and API compatibility unverified.
- Preflight root LICENSE warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-crispr-screens-combinatorial-screens). Original author:
GPTomics; upstream GPTomics/bioSkills path `crispr-screens/combinatorial-screens`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
