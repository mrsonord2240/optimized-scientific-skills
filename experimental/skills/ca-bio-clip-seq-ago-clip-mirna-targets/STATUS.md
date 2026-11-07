# Experimental status: ca-bio-clip-seq-ago-clip-mirna-targets

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-clip-seq-ago-clip-mirna-targets`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Preserved optimized source scripts, methods, provenance and failure references; split direct chimera, targeted UMI and AGO-prediction routes.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Standard AGO peaks do not establish a direct miRNA pairing; seed scan and peak calling are separate workflows.
- Named Hyb database, protocol-specific preprocessing and reference mapping must be prepared independently; exact interfaces require verification.
- TargetScan converter overwrites output paths; use new paths.
- No scientific execution, installs or testing; source tests reviewed but not run.
- Manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-clip-seq-ago-clip-mirna-targets). Original author:
GPTomics; upstream GPTomics/bioSkills path `clip-seq/ago-clip-mirna-targets`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
