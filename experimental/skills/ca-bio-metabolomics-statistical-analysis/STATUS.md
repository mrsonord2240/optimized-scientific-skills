# Experimental status: ca-bio-metabolomics-statistical-analysis

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-metabolomics-statistical-analysis`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Intent routes, file-based Welch wrapper and retained PCA/model examples.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No computations or scientific runs executed.
- Differential wrapper omits features with fewer than three observations per group; CI, paired, covariate and mixed models require adaptation.
- Supervised model demonstration is synthetic and not full external validation.
- No source LICENSE file; repository license evidence required.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-metabolomics-statistical-analysis). Original author:
GPTomics; upstream GPTomics/bioSkills path `metabolomics/statistical-analysis`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
