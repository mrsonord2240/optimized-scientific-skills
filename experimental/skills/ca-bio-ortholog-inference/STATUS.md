# Experimental status: ca-bio-ortholog-inference

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-ortholog-inference`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Intent router, task-specific routes and file-based client entry points; retained source clients and examples.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No live service calls or scientific workflows executed; APIs, schemas, versions and terms require verification.
- Comparison example does not normalize identifiers or automatically intersect calls.
- eggNOG scripted API and de novo inference are unsupported by the source.
- OMA helper conflates 404 and empty results.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-ortholog-inference). Original author:
GPTomics; upstream GPTomics/bioSkills path `database-access/ortholog-inference`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
