# Experimental status: ca-bio-splicing-quantification

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-splicing-quantification`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Routed BAM events, TPM PSI, clusters, LSVs and special events; retained importable provider helper and topical references.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Normalization only: scientific commands unexecuted; no behavioral audit.
- Source MAJIQ V3 command uses apparently older-generation flags and files; marked unverified in route, requires version-specific audit.
- Source failure-modes reference recommends classical SUPPA2 at very low replication, conflicting with source differential-splicing skill; needs scientific review before use.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-splicing-quantification). Original author:
GPTomics; upstream GPTomics/bioSkills path `alternative-splicing/splicing-quantification`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
