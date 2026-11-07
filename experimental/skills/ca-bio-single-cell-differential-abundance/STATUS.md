# Experimental status: ca-bio-single-cell-differential-abundance

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-single-cell-differential-abundance`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

File-based Milo and scCODA workflows with replicate checks and source compositional and design guidance.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific analyses, installations or API checks.
- At least two samples per condition is a floor, not a power guarantee; covariates and pairing require adaptation.
- scCODA automatic reference and HMC convergence require explicit review; output contains credible effects rather than the complete model.
- sccomp and propeller require validated input-specific workflows.
- Skill-root LICENSE warning: cite repository license in manifest.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-single-cell-differential-abundance). Original author:
GPTomics; upstream GPTomics/bioSkills path `single-cell/differential-abundance`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
