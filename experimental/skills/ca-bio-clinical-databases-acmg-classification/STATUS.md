# Experimental status: ca-bio-clinical-databases-acmg-classification

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-clinical-databases-acmg-classification`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Retained source framework, MIT license and helper; added file-based non-patient educational evidence ledger and routes.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific execution, service calls, patient data or installation.
- Criteria and API currentness not verified; evidence requires current authoritative specifications and qualified review.
- Ledger cannot assess disease, transcript, evidence validity or clinical applicability.
- Somatic method is context-dependent, not a turnkey clinical classifier.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-clinical-databases-acmg-classification). Original author:
GPTomics; upstream GPTomics/bioSkills path `clinical-databases/acmg-classification`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
