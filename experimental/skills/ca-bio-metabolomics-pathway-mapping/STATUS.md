# Experimental status: ca-bio-metabolomics-pathway-mapping

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-metabolomics-pathway-mapping`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Task routes, retained source scripts and moved local ORA BH adjustment before output hit filter.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific calculations, network calls or analyses executed.
- Source mummichog sets mpt without RT but v2 requires RT; marked incomplete.
- Local ORA relabels generic pathway IDs as hsa without confirming organism membership; flagged.
- Remote API sends user compound data; disclosure required.
- No source LICENSE file; repository license evidence required.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-metabolomics-pathway-mapping). Original author:
GPTomics; upstream GPTomics/bioSkills path `metabolomics/pathway-mapping`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
