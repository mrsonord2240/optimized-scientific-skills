# Experimental status: ca-bio-causal-genomics-pleiotropy-detection

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-causal-genomics-pleiotropy-detection`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Intent router and separate sensitivity, PRESSO, SIMEX, CAUSE, weak-instrument and input-dependent routes; reference context and original GPTomics/MIT attribution retained.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Normalization only: no scientific workflow, package installation or live API verification was performed.
- LCV, LHC-MR, bidirectional MR, contamination mixture and MRMix need study-specific inputs and adaptation; no turnkey wrapper is claimed.
- Steiger requires sufficient frequency and sample-size fields; CAUSE requires externally selected and pruned SNP lists; source notes CAUSE and loo version incompatibility.
- Isq export and package APIs may vary by installed version.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-causal-genomics-pleiotropy-detection). Original author:
GPTomics; upstream GPTomics/bioSkills path `causal-genomics/pleiotropy-detection`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
