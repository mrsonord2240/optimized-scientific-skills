# Experimental status: ca-bio-atac-seq-single-cell-atac

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-atac-seq-single-cell-atac`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Routed Signac, ecosystems, doublets and multiome; retained source implementation and references with MIT notice.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Normalization only; all scientific workflows unexecuted. Script fixed to hg38 and does not run doublet detection or expose output prefix in CLI.
- ARC blacklist ratio is computed from the peak matrix and may not behave like singlecell.csv blacklist fraction; source warns it may remove no cells.
- ArchR, SnapATAC2, WNN, PEAKVI and AMULET reference workflows need their own tool environments and were not run here.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-atac-seq-single-cell-atac). Original author:
GPTomics; upstream GPTomics/bioSkills path `atac-seq/single-cell-atac`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
