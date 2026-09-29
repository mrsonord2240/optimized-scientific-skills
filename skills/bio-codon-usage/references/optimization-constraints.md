# Max-CAI optimization constraints

## What Biopython optimizes

`CodonAdaptationIndex.optimize()` selects the single highest-weight synonymous
codon for each amino acid. The resulting sequence is a textbook max-CAI draft,
not a multi-objective expression design.

```python
from codon_utils import optimize_dna_table_aware

optimized, source_check, optimized_check = optimize_dna_table_aware(
    cai, query, table_id=2, strict=True
)
assert optimized.translate(table=2) == query.translate(table=2)
```

The wrapper translates and verifies with the same selected table. This avoids
Biopython 1.85's standard-code translation in the direct DNA optimization
path. `strict=True` raises if two synonymous codons tie for the highest weight;
in 1.85, `strict=False` silently chooses one. Protein preservation is necessary
but is not sufficient evidence that a construct will express well.

## Screen the draft before synthesis

- **Translation ramp:** slow codons in roughly the first 30-50 codons can
  affect ribosome spacing. Flattening the ramp may reduce yield or increase
  misfolding.
- **5' RNA secondary structure:** strong folding near the start codon can
  impede initiation, even if CAI rises.
- **GC extremes:** very high GC can stabilize unwanted structures; very low GC
  introduces other stability and synthesis risks.
- **Cryptic regulatory elements:** synonymous swaps can create splice sites,
  Shine-Dalgarno-like sequences, polyadenylation signals, restriction sites,
  or destabilizing elements.
- **Codon-pair effects:** CAI scores codons independently and ignores adjacent
  pairs.

Report these screens, their software and parameters, and any compromise made
against the max-CAI sequence. High CAI does not establish high expression.

## tRNA Adaptation Index

tAI is a supply-side alternative that weights codons by tRNA gene copy number
as a proxy for abundance, adjusted for third-position wobble efficiency. It
may track elongation behavior better where copy number is a useful abundance
proxy. Biopython does not implement tAI; use a separately validated tool such
as the R `tAI` package or a documented implementation.

