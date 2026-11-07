# Experimental status: ca-bio-sra-data

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-sra-data`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Task routes with public access and size boundaries; retained generic source helpers.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No live services, downloads or scientific workflows executed.
- Source Entrez demonstration omitted because it embeds specific queries; generic Entrez discovery is routed to sibling skill.
- ENA batch script can overwrite files and exit zero with failures; route requires empty output and failed.txt review.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-sra-data). Original author:
GPTomics; upstream GPTomics/bioSkills path `database-access/sra-data`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
