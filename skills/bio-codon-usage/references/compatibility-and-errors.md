# Biopython compatibility and common errors

## Biopython 1.82 API migration

`Bio.SeqUtils.CodonUsage` and `Bio.SeqUtils.CodonUsageIndices` were removed in
Biopython 1.82. Use:

```python
from Bio.SeqUtils import CodonAdaptationIndex
```

| Removed interface (through 1.81) | Current interface (1.82+) |
|---|---|
| construct `CodonAdaptationIndex()` then call `generate_index(fasta)` | construct `CodonAdaptationIndex(reference_seqs, table=...)` |
| `cai.cai_for_gene(seq)` | `cai.calculate(seq)` |
| `cai.set_cai_index(weights)` | `cai.update(weights)` because the object is a `dict` subclass |
| bundled `SharpEcoliIndex` | no replacement; supply expression-biased reference CDS |
| `cai.print_index()` | iterate `cai.items()` |

Exact Biopython 1.85 behavior verified for this skill:

- `calculate()` hard-excludes ATG and TGG. A stop codon present in the mapping
  contributes; only a missing TAA, TAG, or TGA is silently skipped.
- Unobserved codons begin with a 0.5 pseudocount, then that value is normalized
  within the synonymous family. For example, unobserved GCC weighs 0.05 when
  GCT is observed ten times.
- `strict=True` rejects an equally preferred synonymous tie. Despite its
  docstring, 1.85 `strict=False` chooses one without emitting a warning.
- Empty or all-hard-excluded queries reach a zero denominator in raw
  `calculate()`. Use `calculate_cai_guarded(...)` for an actionable error and
  an exclusion report.
- Construction and calculation uppercase sequence input internally. Illegal
  or trailing query codons otherwise surface raw `TypeError` exceptions, which
  is why every routed call validates first.
- `optimize(..., seq_type="DNA")` uses standard-code translation even when the
  index has a nonstandard table. Use `optimize_dna_table_aware(...)` instead.

## Troubleshooting

| Symptom | Likely cause | Route |
|---|---|---|
| `ImportError` for `CodonUsage` | removed module | import `CodonAdaptationIndex` directly from `Bio.SeqUtils` |
| missing `generate_index` | obsolete class API | pass reference sequences to the constructor and use `calculate` |
| plausible CAI from a shifted CDS | wrong reading-frame origin | establish codon position 1 and validate the CDS before scoring |
| `ZeroDivisionError` from `calculate()` | query is empty or contains only ATG/TGG | use `calculate_cai_guarded(...)`; revise the query rather than inventing a score |
| CAI near 1 for most genes | uninformative reference set | use highly expressed genes from the exact target host |
| `ValueError` about equally preferred codons | `optimize(strict=True)` encountered a tie | curate the weights or deliberately use `strict=False` |
| alternate-code protein changes after DNA optimization | Biopython's DNA path translated with the standard table | use `optimize_dna_table_aware(...)` and verify with the selected table |
| high CAI but poor measured expression | max-CAI omitted other design constraints | screen the translation ramp, RNA structure, GC, cryptic elements, and codon pairs |

When the installed version differs from the reference version, use
`pip show biopython` and Python `help()` or signature introspection before
changing a call.

