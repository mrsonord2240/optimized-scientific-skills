# Experimental status: ca-bio-crispr-screens-batch-correction

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-crispr-screens-batch-correction`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Split diagnosis, ComBat, control anchoring, RUV, SVA and MLE into intent routes; retained source guidance and extracted complete file-based workflows. GPTomics MIT attribution retained.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- Normalization only; no scientific workflow or installation executed.
- Source guidance inconsistently suggests a covariate fixes complete confounding; candidate flags non-identifiability.
- Diagnostic and ComBat scripts silently select matching samples; verify metadata and sample names.
- Package and API compatibility, numerical performance and downstream analysis remain unverified.
- Preflight warns that root LICENSE is absent; repository MIT license evidence required in publication manifest.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-crispr-screens-batch-correction). Original author:
GPTomics; upstream GPTomics/bioSkills path `crispr-screens/batch-correction`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
