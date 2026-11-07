# Experimental status: ca-bio-chipseq-allele-specific-binding

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-chipseq-allele-specific-binding`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Created intent router and task routes; retained source method references, runners and tests.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No analyses, tests, installations, or external tool invocations were performed.
- AlleleSeq remains external-only; WASP requires assembling the source multi-command workflow and variant tables.
- Source version contracts and API claims not independently verified during normalization.
- Repository-root license evidence is outside candidate.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-chipseq-allele-specific-binding). Original author:
GPTomics; upstream GPTomics/bioSkills path `chip-seq/allele-specific-binding`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
