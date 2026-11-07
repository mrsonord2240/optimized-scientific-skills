# Experimental status: ca-bio-interaction-databases

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-interaction-databases`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

User-file resource query wrapper and task routes around preserved source clients and license.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Commercial OmniPath screening is source-level only; computed fields and per-source terms require independent review. Aggregate defaults to academic scope, source caller identity, undirected edges and no self-interactions.
- IntAct, Reactome, HuRI and HuMAP remain selection guidance; no source runnable client is supplied.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-interaction-databases). Original author:
GPTomics; upstream GPTomics/bioSkills path `database-access/interaction-databases`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
