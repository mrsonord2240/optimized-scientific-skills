# Experimental status: ca-bio-workflows-metabolomics-pipeline

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-workflows-metabolomics-pipeline`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Preserved three source scripts and built file-based extraction, import and QC wrappers, with separate interpretation handoffs.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Annotation, guarded statistics and pathway mapping require external component workflows and adaptation; no complete end-to-end wrapper.
- Blank and D-ratio filters are not implemented by source QC script; metadata beyond sample_class must be reconciled after dropped samples.
- No scientific workflow or network access executed; runtime APIs unverified.
- Manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-workflows-metabolomics-pipeline). Original author:
GPTomics; upstream GPTomics/bioSkills path `workflows/metabolomics-pipeline`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
