# Codon-usage metrics and method details

## Required imports

```python
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Data.CodonTable import standard_dna_table
from Bio.Seq import Seq
from Bio.SeqUtils import CodonAdaptationIndex, GC123
from codon_utils import (
    calculate_cai_guarded,
    count_codons,
    optimize_dna_table_aware,
    validate_cds,
)
```

Import only the components needed for the selected analysis.

## Codon Adaptation Index

CAI measures how closely a query's codon usage resembles an
expression-biased reference set for a target host. Build the index from highly
expressed target-host genes; modern Biopython does not bundle a reference.

```python
reference_seqs = list(SeqIO.parse("highly_expressed_genes.fasta", "fasta"))
cai = CodonAdaptationIndex(reference_seqs, table=standard_dna_table)

query = Seq("ATGAAACGTGCTGAAGCTAAATAA")
score = calculate_cai_guarded(cai, query, table_id=1)
print(f"CAI: {score.score:.3f}; included={score.included_codons}; excluded={score.excluded_codons}")
print(cai["GCT"])
```

The index is a codon-to-relative-adaptiveness mapping and a `dict` subclass.
Calling `cai.update({"GCT": 0.9})` deliberately overrides a weight.

For each synonymous family, CAI's relative-adaptiveness weight is the codon's
RSCU divided by the largest RSCU in the family. The final CAI is the geometric
mean of included weights.

## Relative Synonymous Codon Usage

RSCU divides a codon's observed count by the count expected if the amino
acid's synonymous codons were used uniformly. It normalizes away amino-acid
composition.

```python
def calculate_rscu(seq, table_id=1, policy="strict"):
    table = CodonTable.unambiguous_dna_by_id[table_id]
    counts, validation = count_codons(
        seq, table_id=table_id, policy=policy
    )
    by_amino_acid = {}
    for codon, amino_acid in table.forward_table.items():
        by_amino_acid.setdefault(amino_acid, []).append(codon)

    result = {}
    for codons in by_amino_acid.values():
        total = sum(counts.get(codon, 0) for codon in codons)
        expected = total / len(codons) if codons else 0
        for codon in codons:
            result[codon] = counts.get(codon, 0) / expected if expected else 0
    return result, validation
```

## Effective number of codons is not implemented

Nc is a reference-free codon-bias measure, conventionally ranging from about
20 for complete bias to 61 for no bias. The provider's helper summed
per-amino-acid reciprocal homozygosities, but that quantity is not Wright's
standard class-averaged Nc. It diverged from codonW on bounded comparisons and
has therefore been removed instead of being presented under the Nc name.

Use a separately validated Wright/codonW-compatible implementation when Nc is
required. Record the exact estimator, version, genetic code, and input policy.
Do not compare the removed approximation across studies or convert its values
to Nc by clipping them to 20-61.

## GC at codon positions

```python
gc_total, gc_pos1, gc_pos2, gc_pos3 = GC123(seq)
print(f"GC3: {gc_pos3:.1f}%")
```

`GC123` returns four percentages from 0 to 100. `gc_fraction` uses a 0-to-1
scale. GC3 often reflects genome-wide compositional bias.

## Genetic codes

```python
table = CodonTable.unambiguous_dna_by_id[1]
print(table.start_codons, table.stop_codons)
print(table.forward_table["ATG"])
```

| ID | Name | Typical use |
|---|---|---|
| 1 | Standard | most nuclear genomes |
| 2 | Vertebrate mitochondrial | human and mouse mitochondria |
| 4 | Mold/protozoan mitochondrial | fungal and protozoan mitochondria |
| 5 | Invertebrate mitochondrial | insect and worm mitochondria |
| 11 | Bacterial/plastid | bacteria and chloroplasts |

Pass the matching table to `CodonAdaptationIndex`, `validate_cds`,
`calculate_cai_guarded`, and `optimize_dna_table_aware`; do not assume table 1.
The table-aware optimizer translates to protein with the selected table before
asking Biopython to choose preferred codons, then verifies translation again
under that same table.

## Metric summary

| Metric | Typical range | Reference required | Interpretation |
|---|---:|---|---|
| CAI | 0-1 | highly expressed genes from the target host | larger means closer to that reference |
| RSCU | 0-N | none beyond the analyzed sequence | 1 means uniform synonymous use |
| GC3 | 0-100% | none | GC at codon position 3 |
| tAI | 0-1 | tRNA gene-copy or abundance information | larger means greater inferred tRNA supply |

Nc is intentionally absent from this implemented-metric table; see the
non-implementation boundary above.

