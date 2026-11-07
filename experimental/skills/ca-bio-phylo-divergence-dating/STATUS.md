# Experimental status: ca-bio-phylo-divergence-dating

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-phylo-divergence-dating`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Split calibrated prior/posterior dating, temporal signal and infinite-sites analysis; retained source helpers and GPTomics MIT attribution.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Randomization helper uses fixed filenames and IQ-TREE -redo; use a fresh dedicated directory. It compares point estimates against a randomized range, not CI overlap, and does not check subprocess exit before parsing.
- MCMCTree setup with two inputs immediately runs a full pipeline and uses source default priors; review before use.
- Root-to-tip optimized R-squared is exploratory and does not prove temporal signal.
- No tree inference, simulations, analyses or installations.
- Root LICENSE warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-phylo-divergence-dating). Original author:
GPTomics; upstream GPTomics/bioSkills path `phylogenetics/divergence-dating`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
