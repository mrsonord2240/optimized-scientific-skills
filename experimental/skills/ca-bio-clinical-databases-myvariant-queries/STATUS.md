# Experimental status: ca-bio-clinical-databases-myvariant-queries

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-clinical-databases-myvariant-queries`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Added file-based public annotation and gene search entry points with version sidecar, preserved original helper and historical source context.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No live service calls, workflow execution, installations or patient data.
- Service schemas, index build, source versions and clinical criteria unverified.
- Gene search returns limited hits, and source predicates may match different RCVs; record total and examine conditions.
- Source AlphaMissense statements conflict with later source context; current expert guidance needed.
- Skill-root LICENSE warning: cite repository license in manifest.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-clinical-databases-myvariant-queries). Original author:
GPTomics; upstream GPTomics/bioSkills path `clinical-databases/myvariant-queries`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
