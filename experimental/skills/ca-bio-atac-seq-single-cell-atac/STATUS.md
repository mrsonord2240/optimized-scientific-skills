# Experimental status: ca-bio-atac-seq-single-cell-atac

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft revises `bio-atac-seq-single-cell-atac`. Its earlier audit
score does not certify this changed tree. Do not present the draft as Production
Ready or Marketplace-accepted.

## Structural changes

Routed Signac, ArchR, SnapATAC2, cluster peak calling and WNN; added file-based entry points from source workflows.

The 2026-11-07 route review replaced package-only stubs with task entry points
where the source provides a workflow. Input-dependent or incomplete source
methods are identified explicitly. New wrappers remain unexecuted.

## Validation and remaining work

The revised tree passed offline structural preflight and applicable syntax
checks. Capability coverage and command interfaces were reviewed. Tooling
preparation and meaningful behavioral validation remain before relying on results.

- Normalization only; all scientific workflows and APIs remain unexecuted. ArchR and SnapATAC2 entry points assume hg38 and require installed matching tools and MACS3. SnapATAC2 gene activity remains an input-dependent source reference and is not saved by the CLI.
- Signac script assumes hg38 and does not run doublet detection or expose an output prefix. ARC blacklist ratio comes from the peak matrix and may not behave like the ATAC singlecell.csv fraction.
- Cell-cycle correction, sex QC, PEAKVI mapping, chromBPNet and AMULET remain source-specific workflows requiring additional inputs and adaptation.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-atac-seq-single-cell-atac).
Original author: GPTomics; upstream GPTomics/bioSkills path `atac-seq/single-cell-atac`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This modified derivative retains
the repository [MIT license](../../../LICENSE). The
[manifest](../../CANDIDATES.json) records file hashes, coverage, validation notes
and previous published identities.
