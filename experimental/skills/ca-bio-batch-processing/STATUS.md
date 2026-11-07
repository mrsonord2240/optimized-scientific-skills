# Experimental status: ca-bio-batch-processing

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft was published on 2026-11-07. It revises
`bio-batch-processing`. The original Skill's historical audit score does not certify
this changed tree. Do not present this draft as Production Ready or
Marketplace-accepted.

## Structural changes

Intent routes, source recipes and reader references, reusable CLI over source functions; MIT GPTomics attribution retained.

## Validation and remaining work

The normalized tree passed the project's offline structural preflight. Syntax
checks and capability-to-route coverage were reviewed; scientific workflows
were not executed. Tooling preparation and a meaningful behavioral audit are
the next phases before relying on results.

- No scientific runs, library imports, regression execution or dependency/version verification.
- Plain-text core paths do not transparently handle gzip; use pyfastx or adapt reader.
- Split by count, conversion, merge and summary may overwrite outputs; select fresh paths.
- Retained source recipes require task-specific edits, including source tracking and multi-file index.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-batch-processing). Original author:
GPTomics; upstream GPTomics/bioSkills path `sequence-io/batch-processing`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This is a modified derivative
under the repository [MIT license](../../../LICENSE). See the
[candidate manifest](../../CANDIDATES.json) for source identity, reviewed
capabilities, validation notes, and the exact normalized file hashes.
