---
name: bio-codon-usage
description: Analyze codon usage and calculate CAI (Codon Adaptation Index) and RSCU with Biopython, and produce table-aware naive max-CAI codon-optimized sequences. Use when scoring a gene's codon bias against a host, optimizing a CDS for heterologous expression, or studying synonymous codon selection.
tool_type: python
primary_tool: Bio.SeqUtils.CodonAdaptationIndex
license: MIT
author: GPTomics
---

# Codon Usage

Analyze in-frame coding sequences, quantify synonymous-codon bias, score a gene
against an expression-biased host reference, or generate a max-CAI draft for
further design screening.

## Version compatibility

The routed examples target Biopython 1.83+. Before reusing them, inspect the
installed version with `pip show biopython` and use Python `help(...)` for any
changed signature. If an import, attribute, or type error occurs, inspect the
installed API rather than retrying an obsolete call.

Modern Biopython no longer provides `Bio.SeqUtils.CodonUsage` or
`CodonUsageIndices`. Import `CodonAdaptationIndex` directly from
`Bio.SeqUtils`. See [compatibility and errors](references/compatibility-and-errors.md)
for the 1.82 migration table and failure routes.

## Choose the task

- For codon counts, frequencies, amino-acid usage, and GC at each codon
  position, run or adapt [basic_analysis.py](scripts/basic_analysis.py).
- For RSCU and low-RSCU codons, run or adapt
  [rscu_analysis.py](scripts/rscu_analysis.py).
- For host-relative CAI and a naive max-CAI draft, start with
  [cai_optimization.py](scripts/cai_optimization.py).
- For alternate genetic codes, metric definitions, and implementation patterns,
  read [metrics and methods](references/metrics-and-methods.md).
- Before treating an optimized sequence as design-ready, follow
  [optimization constraints](references/optimization-constraints.md).

## Validate every input before analysis

1. Run `validate_cds(...)` from
   [codon_utils.py](scripts/codon_utils.py) before every routed analysis.
2. Use `policy="strict"` for scientific results. It rejects empty input,
   partial triplets, ambiguity, a table-invalid start or terminal stop, and
   internal stops. Use `policy="permissive"` only for explicit cleanup: it
   discards whole ambiguous codons and a trailing partial codon and records
   every discard in `CDSValidationResult.discarded`.
3. Confirm the sequence is a CDS in the intended orientation and that base 1
   is codon position 1. The start check detects many shifted inputs but cannot
   establish orientation or biological annotation by itself.
4. Select the correct genetic code; do not assume the standard table for
   mitochondrial or other nonstandard genes.
5. For CAI, supply CDS from highly expressed genes of the target organism,
   such as ribosomal proteins and elongation factors. Do not substitute a
   whole-genome average or a different organism.

Do not call Biopython metrics on raw input. The shared validator makes length,
symbols, and selected-table boundaries executable, while orientation and the
annotated reading-frame origin remain caller responsibilities.

## Run the requested analysis

### Count codons or calculate RSCU

Count non-overlapping triplets from the established reading frame. For RSCU,
group sense codons by encoded amino acid, then divide each observed count by
the within-family count expected under uniform synonymous usage. Interpret
RSCU 1 as unbiased within the family, values above 1 as over-used, and values
below 1 as under-used.

Use the routed helpers as executable starting points. Report the genetic code,
sequence length, analyzed codon count, and any discarded or invalid symbols.

### Build and apply a CAI index

Validate every reference and query, construct
`CodonAdaptationIndex(reference_seqs, table=...)` from an expression-biased
reference set for the target host, then call `calculate_cai_guarded(...)`.
The wrapper checks table consistency, prevents the raw zero-denominator error,
and reports codons excluded by Biopython 1.85. Report the reference-set
definition and genetic code with the score; CAI without that context has no
biological meaning.

Treat `CodonAdaptationIndex` as the codon-to-relative-adaptiveness mapping. It
is a `dict` subclass and may be inspected or deliberately updated. Do not
claim that modern Biopython bundles a Sharp E. coli index.

### Generate a max-CAI draft

Call `optimize_dna_table_aware(...)`, which validates the DNA, translates it
with the selected genetic-code table, optimizes the protein sequence, then
retranslates and verifies it with the same table. Do not call Biopython 1.85
`optimize(..., seq_type='DNA')` for a nonstandard table: that path translates
DNA with the standard code. `strict=True` raises when synonymous codons tie for
highest weight; in Biopython 1.85, `strict=False` silently chooses one.

Never present this output as an expression-optimized final construct. It
maximizes a single-codon frequency objective while ignoring translation-ramp
shape, 5' RNA structure, GC extremes, cryptic regulatory elements, and codon
pairs. Screen those constraints before synthesis.

## Report with bounded interpretation

- CAI is host-reference-relative; higher means closer to the supplied
  expression-biased reference, not guaranteed higher expression.
- RSCU removes amino-acid-composition effects within synonymous families.
- GC123 returns percentages from 0 to 100, unlike `gc_fraction`.
- This skill does not calculate Nc. The provider's per-amino-acid reciprocal-
  homozygosity sum was removed because it is not Wright/codonW Nc and can
  diverge materially. Use a separately validated standard implementation for
  Nc, and record the exact estimator and version; never compare the removed
  approximation across studies.
- tAI is a supply-side alternative based on tRNA gene copy number and wobble
  efficiency; it is not implemented in Biopython.

State assumptions, reference provenance, frame and table checks, software
version, and any limitations alongside results.

## Resources

- [Metrics and methods](references/metrics-and-methods.md) — CAI, RSCU,
  GC123, genetic-code tables, the Nc exclusion, and implementation details.
- [Optimization constraints](references/optimization-constraints.md) —
  translation ramp, RNA structure, GC, cryptic elements, codon pairs, tAI,
  and result-screening requirements.
- [Compatibility and errors](references/compatibility-and-errors.md) —
  Biopython API migration and symptom-driven troubleshooting.
- [Provenance](references/provenance.md) — upstream identity, preserved script
  blobs, license, and scientific references.

