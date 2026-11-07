# Experimental status: ca-bio-crispr-screens-jacks-analysis

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-crispr-screens-jacks-analysis`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Split joint fit, efficacy summary, and reference reuse; retained source helpers, conditional guidance and GPTomics MIT attribution.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Bundled wrapper does not expose reffile; use source JACKS CLI.
- Source internal prior/convergence overrides rely on non-public version-specific APIs and are not verified.
- Result example uses fixed filenames and first cell line; CLI and dependencies not executed.
- Root LICENSE preflight warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-crispr-screens-jacks-analysis). Original author:
GPTomics; upstream GPTomics/bioSkills path `crispr-screens/jacks-analysis`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
