# Experimental status: ca-bio-batch-downloads

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-batch-downloads`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Task routes with user-input history wrapper, source retry worker and known-ID guidance.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Normalization only; no network or scientific workflow execution. History offset resume assumes result ordering does not change; old checkpoints lack a strong boundary.
- Source EPost limits and WebEnv lifetime are historical assertions pending current primary documentation verification; do not treat as current policy.
- Source known-ID example executes a demo on import; reference explicitly warns.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-batch-downloads). Original author:
GPTomics; upstream GPTomics/bioSkills path `database-access/batch-downloads`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
