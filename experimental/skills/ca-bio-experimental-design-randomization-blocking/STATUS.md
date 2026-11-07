# Experimental status: ca-bio-experimental-design-randomization-blocking

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-experimental-design-randomization-blocking`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

File-based blocked assignment and experimental-unit summary routes, with source design and model guidance.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific runs, simulations, package installs or API checks.
- Two-label blocked assignment is not a general constrained or sequential randomizer; odd blocks favor first label by one.
- Unit aggregation is not universally appropriate; source mixed-model designs require study-specific adaptation.
- Skill-root LICENSE warning: cite repository license in manifest.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-experimental-design-randomization-blocking). Original author:
GPTomics; upstream GPTomics/bioSkills path `experimental-design/randomization-blocking`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
