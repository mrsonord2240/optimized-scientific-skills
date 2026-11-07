# Experimental status: ca-bio-vcf-manipulation

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-vcf-manipulation`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Intent routes and focused references preserve source bcftools operations, missingness, sample, build, allele and index caveats; GPTomics MIT credit retained.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Source reference assumes bcftools 1.19+ and cyvcf2 0.30+; installed APIs and scientific workflows not executed.
- Position-only overlap ignores REF/ALT. Structural-variant consensus and gVCF joint genotyping remain upstream-only, not implemented.
- Source helper may create indexes beside inputs and requires a fresh output directory.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-vcf-manipulation). Original author:
GPTomics; upstream GPTomics/bioSkills path `variant-calling/vcf-manipulation`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
