# Experimental status: ca-bio-proteomics-protein-inference

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-proteomics-protein-inference`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Split Basic, EPIFANY, picked-group and external CLI paths; source helpers and GPTomics MIT attribution retained.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Picked table helper uses exact group-set pairing and fixed 1% cutoff; mismatched memberships remain unpaired.
- EPIFANY greedy resolution does not necessarily remove subsumable groups.
- External output and package APIs not executed.
- No scientific execution or installations.
- Root LICENSE warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-proteomics-protein-inference). Original author:
GPTomics; upstream GPTomics/bioSkills path `proteomics/protein-inference`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
