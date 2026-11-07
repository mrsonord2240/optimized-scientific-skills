# Experimental status: ca-bio-single-cell-scatac-analysis

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-single-cell-scatac-analysis`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

File-based prepared-object LSI and retained chromVAR workflow; source QC, peaks and alternate framework context routed.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific analyses, installations or API checks.
- Prepared-object LSI examines only component 1 automatically; inspect all components and tailor QC.
- chromVAR source helper densifies counts and requires substantial memory; cell-level differential motif results are descriptive without replicate-aware design.
- ArchR, SnapATAC2, peak calling and doublet methods require adapted workflows.
- Skill-root LICENSE warning: cite repository license in manifest.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-single-cell-scatac-analysis). Original author:
GPTomics; upstream GPTomics/bioSkills path `single-cell/scatac-analysis`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
