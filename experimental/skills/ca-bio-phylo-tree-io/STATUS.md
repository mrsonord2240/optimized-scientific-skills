# Experimental status: ca-bio-phylo-tree-io

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-phylo-tree-io`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Split original tree parsing, annotation extraction, format conversion and background into focused routes and references; adapted synthetic examples to file-based scripts. Conversion refuses existing output and identical input/output.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- DendroPy and Bio.Phylo runtime APIs and full annotation round-trip were not executed or verified.
- JSON fallback may stringify metadata values; compare the original.
- Plain conversion may lose annotations and requires a new destination.
- No analyses, installs, simulations or inference were run.
- Root LICENSE warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-phylo-tree-io). Original author:
GPTomics; upstream GPTomics/bioSkills path `phylogenetics/tree-io`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
