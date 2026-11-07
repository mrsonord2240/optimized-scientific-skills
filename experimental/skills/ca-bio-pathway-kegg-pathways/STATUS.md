# Experimental status: ca-bio-pathway-kegg-pathways

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-pathway-kegg-pathways`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Intent routes, explicit-universe KEGG ORA wrapper and retained source examples.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No enrichment, network or scientific workflows executed.
- Wrapper uses live KEGG and does not snapshot or report post-mapping coverage; inspect independently.
- GSEA, modules, topology, map overlay and snapshots require adaptation.
- No source LICENSE file; repository license evidence required.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-pathway-kegg-pathways). Original author:
GPTomics; upstream GPTomics/bioSkills path `pathway-analysis/kegg-pathways`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
