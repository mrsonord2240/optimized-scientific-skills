# Experimental status: ca-bio-entrez-link

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-entrez-link`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Routed direct discovery, large union and chained links with user-input wrappers, retaining source helper.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Normalization only; no live data calls. Source linknames and rate limits require current verification.
- History wrapper records source count and history identifiers, but linked count requires a follow-up ESearch and source IDs must be validated separately.
- Source checker uses ESummary schemas and may fail on schema changes; label confirmation cannot prove biological identity.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-entrez-link). Original author:
GPTomics; upstream GPTomics/bioSkills path `database-access/entrez-link`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
