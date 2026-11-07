# Experimental status: ca-bio-clinical-databases-gnomad-frequencies

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-clinical-databases-gnomad-frequencies`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Public-example variant and gene CLI with source helper and historical references; undefined FAF95 requires review and contradictory LOEUF decile is not auto-assigned.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No live service calls, workflow execution, installs or patient data.
- Source release, GraphQL schema and clinical criteria not verified.
- Source presents conflicting v4 LOEUF decile thresholds; automatic decile flag disabled.
- Missing FAF95 may mean no observation, undefined value or build mismatch; qualified review required.
- Bulk, SV/CNV and mtDNA require validated data and workflows.
- Skill-root LICENSE warning: cite repository license in manifest.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-clinical-databases-gnomad-frequencies). Original author:
GPTomics; upstream GPTomics/bioSkills path `clinical-databases/gnomad-frequencies`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
