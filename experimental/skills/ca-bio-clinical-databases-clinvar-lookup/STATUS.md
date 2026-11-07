# Experimental status: ca-bio-clinical-databases-clinvar-lookup

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-clinical-databases-clinvar-lookup`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Routes for local VCF, public summary and allele registry; retained source helper and historical reference context. Unknown review statuses now return null rather than silently mapping to zero stars.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No live calls, workflow execution, package installation, or patient data.
- Source helper imports requests, pandas and cyvcf2 for every mode; installed dependencies and service/schema contracts not verified.
- Condition-specific RCV XML parsing and SCV evidence require a separate validated parser.
- Local VCF source helper returns variant-level fields and may return no result.
- Skill-root LICENSE warning: manifest must cite repository license.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-clinical-databases-clinvar-lookup). Original author:
GPTomics; upstream GPTomics/bioSkills path `clinical-databases/clinvar-lookup`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
