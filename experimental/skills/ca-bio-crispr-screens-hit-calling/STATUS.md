# Experimental status: ca-bio-crispr-screens-hit-calling

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-crispr-screens-hit-calling`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Split consensus, guide support, custom z score and method selection routes; retained source scripts and GPTomics MIT attribution.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Consensus missing values count as noncalls; overlap enrichment is heuristic and cannot establish comparison identity.
- Source statements implying MAGeCK, BAGEL2 and drugZ answer the same question require explicit comparison validation.
- Source inconsistent about whether Chronos training requires CN and amount of validation for tier 1; retain caveats and do not claim consensus validates hits.
- Scientific execution and dependencies unverified.
- Root LICENSE preflight warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-crispr-screens-hit-calling). Original author:
GPTomics; upstream GPTomics/bioSkills path `crispr-screens/hit-calling`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
