# Experimental status: ca-bio-biomart-queries

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-biomart-queries`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Routed ID, coordinate and ortholog queries with a parameterized wrapper around source helper.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Normalization only; no network or scientific execution. Python CLI targets current host; release pinning is the R reference. Schema, endpoint and attribute limits are source-tested historical claims pending live verification.
- The source helper uses pandas TSV parsing and cannot guarantee completeness of a syntactically valid truncated response; validate output.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-biomart-queries). Original author:
GPTomics; upstream GPTomics/bioSkills path `database-access/biomart-queries`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
