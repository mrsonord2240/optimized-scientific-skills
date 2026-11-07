# Experimental status: ca-bio-experimental-design-power-analysis

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-experimental-design-power-analysis`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Added file-based RNASeqPower and proteomics calculations, retained source synthetic bulk and pseudobulk simulations, and extracted design context.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific execution, simulations, package installation or API verification.
- Single-feature approximations do not estimate family-level FDR power.
- PROPER and pseudobulk examples use fixed illustrative priors and need study-specific adaptation; realized FDR estimates have Monte Carlo uncertainty.
- Skill-root LICENSE warning: cite repository license in manifest.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-experimental-design-power-analysis). Original author:
GPTomics; upstream GPTomics/bioSkills path `experimental-design/power-analysis`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
