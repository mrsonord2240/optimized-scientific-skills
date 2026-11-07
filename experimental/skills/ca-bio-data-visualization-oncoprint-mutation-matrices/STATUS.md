# Experimental status: ca-bio-data-visualization-oncoprint-mutation-matrices

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-data-visualization-oncoprint-mutation-matrices`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Task routes, retained source executable scripts and cohort interpretation rules.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No plotting or statistical workflows executed.
- R example always runs interactions and raw MAF-count fallback is not normalized TMB; flagged.
- Additional calls and table export require adaptation; source tests reviewed but not run or copied.
- No source LICENSE file; repository license evidence required.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-data-visualization-oncoprint-mutation-matrices). Original author:
GPTomics; upstream GPTomics/bioSkills path `data-visualization/oncoprint-mutation-matrices`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
