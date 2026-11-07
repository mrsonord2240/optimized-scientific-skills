# Experimental status: ca-bio-pose-validation

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-pose-validation`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Routed pose checks, strain, batch and planarity; retained source functions and removed fixture-only demo entry point.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Normalization only; no scientific execution or behavioral audit.
- PB-valid aggregation may be misleading if expected boolean checks are absent; route requires checking configuration and loading status.
- Reference conformer optimization convergence and reproducibility need behavioral review.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-pose-validation). Original author:
GPTomics; upstream GPTomics/bioSkills path `chemoinformatics/pose-validation`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
