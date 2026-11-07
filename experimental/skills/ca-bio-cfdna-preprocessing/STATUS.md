# Experimental status: ca-bio-cfdna-preprocessing

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft revises `bio-cfdna-preprocessing`. Its earlier audit
score does not certify this changed tree. Do not present the draft as Production
Ready or Marketplace-accepted.

## Structural changes

Intent router, file-based CLI wrapper, focused routes and retained source references, implementation and tests.

The final 2026-11-07 route review clarified that the QC writer overwrites an existing JSON destination and requires a fresh output path. Script logic was unchanged; no scientific workflows were run.

## Validation and remaining work

The revised tree passed offline structural preflight and applicable syntax
checks. Capability coverage and command interfaces were reviewed. Tooling
preparation and meaningful behavioral validation remain before relying on results.

- No science, external tools or tests executed; installed tool version and flag semantics require verification.
- Minimal sWGS and native-end fragmentomics paths are guidance, not complete pipelines.
- Fixed source wrapper thresholds may need assay-specific adaptation; intermediate filenames can overwrite existing files.
- Repository-root license evidence is outside candidate.
- QC JSON output is overwritten if it exists; route requires a fresh output path.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-cfdna-preprocessing).
Original author: GPTomics; upstream GPTomics/bioSkills path `liquid-biopsy/cfdna-preprocessing`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This modified derivative retains
the repository [MIT license](../../../LICENSE). The
[manifest](../../CANDIDATES.json) records file hashes, coverage, validation notes
and previous published identities.
