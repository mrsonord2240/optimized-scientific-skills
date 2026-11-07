# Experimental status: ca-bio-crispr-screens-screen-qc

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-crispr-screens-screen-qc`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Split count representation, replicate QC, essentialome recovery, CN bias and conditional analyses; retained source helpers and GPTomics MIT attribution.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Representation helper lacks full count validation; fixed example requires editing sample maps.
- Composite source averages unscaled log-depth alongside 0-1 metrics despite describing rescaling; not a validated gate.
- Source presents PR-AUC thresholds as universal, but reference set, class prevalence and context affect interpretation.
- CN focal test can be underpowered and absent significance does not establish absence of bias.
- No scientific execution or installations.
- Root LICENSE preflight warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-crispr-screens-screen-qc). Original author:
GPTomics; upstream GPTomics/bioSkills path `crispr-screens/screen-qc`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
