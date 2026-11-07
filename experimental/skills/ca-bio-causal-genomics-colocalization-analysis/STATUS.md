# Experimental status: ca-bio-causal-genomics-colocalization-analysis

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-causal-genomics-colocalization-analysis`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Normalized the source into progressively routed instructions and reusable helpers.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific runs.
- Harmonized output requires manual conversion into two tables; ABF wrapper does not pass MAF to coloc when sdY omitted, so p-value/MAF or estimated sdY needs separate implementation.
- SharePro lacks a runnable source workflow.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-causal-genomics-colocalization-analysis). Original author:
GPTomics; upstream GPTomics/bioSkills path `causal-genomics/colocalization-analysis`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
