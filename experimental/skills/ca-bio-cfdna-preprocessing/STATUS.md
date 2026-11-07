# Experimental status: ca-bio-cfdna-preprocessing

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-cfdna-preprocessing`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Intent router, file-based CLI wrapper, focused routes and retained source references, implementation and tests.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No science, external tools or tests executed; installed tool version and flag semantics require verification.
- Minimal sWGS and native-end fragmentomics paths are guidance, not complete pipelines.
- Fixed source wrapper thresholds may need assay-specific adaptation; intermediate filenames can overwrite existing files.
- Repository-root license evidence is outside candidate.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-cfdna-preprocessing). Original author:
GPTomics; upstream GPTomics/bioSkills path `liquid-biopsy/cfdna-preprocessing`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
