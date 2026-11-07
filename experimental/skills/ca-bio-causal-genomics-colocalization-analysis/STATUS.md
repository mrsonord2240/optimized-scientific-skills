# Experimental status: ca-bio-causal-genomics-colocalization-analysis

**Normalization only; scientific execution and behavioral audit are pending.**

This partially complete draft revises `bio-causal-genomics-colocalization-analysis`. Its earlier audit
score does not certify this changed tree. Do not present the draft as Production
Ready or Marketplace-accepted.

## Structural changes

Replaced package-only route stubs with file-based task entry points; source-incomplete paths remain explicitly identified.

The 2026-11-07 route review replaced package-only stubs with task entry points
where the source provides a workflow. Input-dependent or incomplete source
methods are identified explicitly. New wrappers remain unexecuted.

## Validation and remaining work

The revised tree passed offline structural preflight and applicable syntax
checks. Capability coverage and command interfaces were reviewed. Tooling
preparation and meaningful behavioral validation remain before relying on results.

- No scientific runs.
- Harmonized output requires manual conversion into two tables; ABF wrapper does not pass MAF to coloc when sdY omitted, so p-value/MAF or estimated sdY needs separate implementation.
- SharePro lacks a runnable source workflow.
- HyPrColoc wrapper and region check are syntax-only; tool execution not verified.
- Region check covers a single supplied position, so the entire locus window requires independent review.
- moloc remains an input-dependent source example requiring adaptation.

## Source and credit

Derived from [the published source Skill](https://github.com/mrsonord2240/optimized-scientific-skills/tree/7346d93559e1e5222dfd462ce4d0c2dac70ef35f/skills/bio-causal-genomics-colocalization-analysis).
Original author: GPTomics; upstream GPTomics/bioSkills path `causal-genomics/colocalization-analysis`
at `d91ed3d563019e649dc854c56ccd62551359488a`. This modified derivative retains
the repository [MIT license](../../../LICENSE). The
[manifest](../../CANDIDATES.json) records file hashes, coverage, validation notes
and previous published identities.
