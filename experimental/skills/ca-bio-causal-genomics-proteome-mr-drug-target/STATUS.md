# Experimental status: ca-bio-causal-genomics-proteome-mr-drug-target

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-causal-genomics-proteome-mr-drug-target`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Intent routes, file-based harmonized cis-MR/PAV panel and source curated phenome script; source reference and decision thresholds retained.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific executions, package installations or API verification.
- Complete cis selection/clumping and coloc require study-specific preparation; correlated-IV route explicitly requires LD sign alignment and adaptation.
- PheWAS makes network calls when run, may fail or skip outcomes, and needs OpenGWAS credentials.
- Skill-root LICENSE warning: publication manifest must cite repository license.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-causal-genomics-proteome-mr-drug-target). Original author:
GPTomics; upstream GPTomics/bioSkills path `causal-genomics/proteome-mr-drug-target`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
