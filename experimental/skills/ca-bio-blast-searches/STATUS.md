# Experimental status: ca-bio-blast-searches

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-blast-searches`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Separated search, RID lifecycle and result parsing; replaced hardcoded network examples with user-file parser.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Normalization only; no remote searches or scientific workflows run. Remote RID retention and service parameter claims are historical source statements pending current verification.
- Source coverage calculation uses HSP alignment length over query length and can be misleading for gaps or translated searches; route calls for coordinate and unit review.
- Bit scores are normalized for scoring scheme but comparisons remain subject to query and scoring context; source reference statements need scientific audit.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-blast-searches). Original author:
GPTomics; upstream GPTomics/bioSkills path `database-access/blast-searches`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
