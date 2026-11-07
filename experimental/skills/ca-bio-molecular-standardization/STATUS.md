# Experimental status: ca-bio-molecular-standardization

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-molecular-standardization`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Separated ChEMBL, explicit RDKit and QSAR workflows; preserved source helper and default choices.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Normalization only; no scientific execution or behavioral audit.
- keep_isotopes_col controls isotope stripping in RDKit branch; ChEMBL branch does not use it to modify processing despite reporting the flag; review mixed-isotope behavior.
- QSAR helper returns aggregate exclusion counts rather than source-row-level statuses.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-molecular-standardization). Original author:
GPTomics; upstream GPTomics/bioSkills path `chemoinformatics/molecular-standardization`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
