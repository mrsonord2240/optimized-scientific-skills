# Experimental status: ca-bio-phylo-species-trees

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-phylo-species-trees`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Split ASTER, multicopy, concordance and conditional SVDQuartets/BPP routes; retained source details and GPTomics MIT attribution.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Source treats quartet symmetry/asymmetry as uniquely identifying ILS or introgression; discordance alone is not conclusive.
- BPP/StarBEAST2 require bespoke inputs; external tool interfaces unverified.
- No inference, simulation, analysis or installation.
- Root LICENSE warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-phylo-species-trees). Original author:
GPTomics; upstream GPTomics/bioSkills path `phylogenetics/species-trees`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
