# Experimental status: ca-bio-workflows-scrnaseq-pipeline

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-workflows-scrnaseq-pipeline`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Adapted source Scanpy and Seurat examples to explicit file inputs and new output directories; split single-sample routes from multi-sample handoffs.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Examples start with filtered data: ambient correction and demultiplexing are separate, multi-sample integration and inference require component workflows.
- Outputs are clustered, not validated cell annotations.
- Flat QC cutoffs and gene namespace are illustrative and require review.
- Scanpy output counts layer is subset to HVGs; retain full counts for pseudobulk.
- Runtime API and scientific behavior unexecuted.
- Manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-workflows-scrnaseq-pipeline). Original author:
GPTomics; upstream GPTomics/bioSkills path `workflows/scrnaseq-pipeline`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
