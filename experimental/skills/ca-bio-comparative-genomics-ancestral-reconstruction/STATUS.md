# Experimental status: ca-bio-comparative-genomics-ancestral-reconstruction

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-comparative-genomics-ancestral-reconstruction`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Separate sequence, discrete and continuous intent routes with model and uncertainty references; retain source examples and add continuous file wrapper. GPTomics MIT license preserved.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific execution or API verification; dependency support claims are source historical statements.
- Discrete script currently prints metadata rather than full node and map results; adapt output persistence.
- Codeml helper functions and Bayesian/gene-content methods need task-specific workflows.
- The original construct-design examples are outside this assignment and were not routed as a task.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-comparative-genomics-ancestral-reconstruction). Original author:
GPTomics; upstream GPTomics/bioSkills path `comparative-genomics/ancestral-reconstruction`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
