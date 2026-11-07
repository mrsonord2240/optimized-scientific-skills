# Experimental status: ca-bio-phylo-tree-manipulation

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-phylo-tree-manipulation`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Separated method and caveat references from routes; implemented file-based outgroup and midpoint rooting, pruning, support collapse, ladderizing and MRCA/distance inspection.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- MAD, MinVar and likelihood commands require installed external tools; arbitrary polytomy resolution requires source-described R or ETE adaptation.
- Rerooting can misassign support; verify bipartitions.
- Newick I/O can lose annotations; use bio-phylo-tree-io for annotation-critical work.
- No workflow executed; runtime and scientific behavior unverified.
- Publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-phylo-tree-manipulation). Original author:
GPTomics; upstream GPTomics/bioSkills path `phylogenetics/tree-manipulation`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
