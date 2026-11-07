# Experimental status: ca-bio-proteomics-peptide-identification

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-proteomics-peptide-identification`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Split CLI search/rescoring, pyOpenMS and concatenated/separate FDR routes; retained source scripts and GPTomics MIT attribution.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Source shell reuses any existing target-decoy FASTA without provenance check and writes headerless filtered lists; fresh output and original Percolator files required.
- Separate-search script prints results but does not export a filtered table.
- FDR requires appropriate decoy null, unique spectrum identity and sufficient data.
- No scientific execution, installations or version checks.
- Root LICENSE warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-proteomics-peptide-identification). Original author:
GPTomics; upstream GPTomics/bioSkills path `proteomics/peptide-identification`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
