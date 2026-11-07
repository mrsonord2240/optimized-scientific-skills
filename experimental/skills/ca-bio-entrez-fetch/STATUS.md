# Experimental status: ca-bio-entrez-fetch

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-entrez-fetch`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Routed full records, summaries, variant and history retrieval; generic user-input file writer wraps Entrez.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Normalization only; no live data calls. Database format combinations and rate constraints are source historical guidance pending current verification.
- Generic wrapper rejects empty and HTML responses but may save an XML error document; verify response and identifiers. Source history worker overwrites its output and lacks checkpointing.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-entrez-fetch). Original author:
GPTomics; upstream GPTomics/bioSkills path `database-access/entrez-fetch`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
