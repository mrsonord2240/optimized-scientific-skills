# Experimental status: ca-bio-experimental-design-sample-size

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-experimental-design-sample-size`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

File-based single-parameter sizing and retained pilot and donor scripts, adding donor alignment validation and explicit PROPER summary printing.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific runs, simulations, installations or package verification.
- PROPER script reports actual FDR but does not enforce its target; review error rates and Monte Carlo uncertainty before selecting a design.
- NaN FDR requires inspecting discovery counts; do not assume a cause.
- Study-specific pilot quality, dispersion, experimental units and failure rates require review.
- Skill-root LICENSE warning: cite repository license in manifest.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-experimental-design-sample-size). Original author:
GPTomics; upstream GPTomics/bioSkills path `experimental-design/sample-size`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
