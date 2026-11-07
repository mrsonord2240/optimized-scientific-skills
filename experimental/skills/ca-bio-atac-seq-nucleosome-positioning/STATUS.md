# Experimental status: ca-bio-atac-seq-nucleosome-positioning

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-atac-seq-nucleosome-positioning`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Routes for diagnostics, occupancy and differential analysis; retained source scripts and MIT notice.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Normalization only; all scientific workflows unexecuted. NucleoATAC accumulator defect and unmaintained Python 2 dependency are source-reported; calls require verification.
- Source slop example lacks strand-aware flag; route explicitly requires strand-aware downstream assignment. scPrinter is named but no tested workflow supplied.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-atac-seq-nucleosome-positioning). Original author:
GPTomics; upstream GPTomics/bioSkills path `atac-seq/nucleosome-positioning`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
