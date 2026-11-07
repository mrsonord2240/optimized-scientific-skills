# Experimental status: ca-bio-proteomics-data-import

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-proteomics-data-import`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Split spectra, MaxQuant and DIA-NN imports; extracted real-file MaxQuant CLI and added DIA-NN TSV export; preserve GPTomics MIT attribution.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- TMT path is source snippet rather than bundled file CLI.
- Missingness correlation alone does not establish MNAR.
- Missing flag columns and duplicate DIA-NN group/run values require review.
- No scientific execution or package version verification.
- Root LICENSE warning: publication manifest must cite repository license evidence.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-proteomics-data-import). Original author:
GPTomics; upstream GPTomics/bioSkills path `proteomics/data-import`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
