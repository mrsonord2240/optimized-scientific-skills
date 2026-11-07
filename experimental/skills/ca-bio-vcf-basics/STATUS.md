# Experimental status: ca-bio-vcf-basics

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-vcf-basics`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Intent routes for viewing, querying, conversion and interpretation; source field semantics and cyvcf2 examples extracted into focused references. GPTomics MIT credit preserved.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific execution or API/version verification; source examples assume cyvcf2 0.30+ and bcftools 1.19+.
- Summary cannot distinguish PASS from unfiltered; use raw record.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-vcf-basics). Original author:
GPTomics; upstream GPTomics/bioSkills path `variant-calling/vcf-basics`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
