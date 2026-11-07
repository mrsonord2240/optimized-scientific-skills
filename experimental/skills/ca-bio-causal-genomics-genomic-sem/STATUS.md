# Experimental status: ca-bio-causal-genomics-genomic-sem

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-causal-genomics-genomic-sem`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Normalized the source into progressively routed instructions and reusable helpers.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific runs.
- Factor example also performs a DWLS-only SNP classification; the separate Q_SNP wrapper supplies ML cross-check, as route specifies.
- Several model routes require adapting source snippets.
- Input-dependent HDL or custom model APIs need runtime verification; no scientific execution.
- Custom and stratified source templates have no complete file-based wrapper and are explicitly described as requiring adaptation.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-causal-genomics-genomic-sem). Original author:
GPTomics; upstream GPTomics/bioSkills path `causal-genomics/genomic-sem`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
