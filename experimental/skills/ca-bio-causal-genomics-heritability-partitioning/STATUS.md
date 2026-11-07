# Experimental status: ca-bio-causal-genomics-heritability-partitioning

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-causal-genomics-heritability-partitioning`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Normalized the source into progressively routed instructions and reusable helpers.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific workflows run; HESS aggregate stage and BOLT source were previously unverified.
- Source CBIIT h2-cts requires documented patch; GCTA PCGC option needs installed-interface verification.
- Some inherited shell paths are unquoted and may not handle spaces.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-causal-genomics-heritability-partitioning). Original author:
GPTomics; upstream GPTomics/bioSkills path `causal-genomics/heritability-partitioning`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
