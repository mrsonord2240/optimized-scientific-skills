# Experimental status: ca-bio-crispr-screens-crispresso-editing

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-crispr-screens-crispresso-editing`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Split single, base/prime, multiple-mode and result workflows; preserved scripts, reference guidance and GPTomics MIT credit.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Source gives inconsistent base-editing window center suggestions; verify coordinates and installed-version semantics.
- Source editing thresholds are conventions and require context-specific interpretation.
- Parser divides by input reads, which must be nonzero; example shell contains fixed sequences.
- CRISPResso APIs, version, outputs and scientific workflows not executed or verified.
- Root LICENSE warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-crispr-screens-crispresso-editing). Original author:
GPTomics; upstream GPTomics/bioSkills path `crispr-screens/crispresso-editing`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
