# Substitution counts

Read when tabulating substitution pairs (or Ti/Tv) from an alignment, or when a bit score depends on the BLOSUM62 lambda. Assumes the normalised `alignment` and the imports from "Required Import and Normalisation" in `SKILL.md`.

## Substitution Counts

**Goal:** Tabulate observed substitution counts from the alignment for evolutionary analysis.

**Approach:** Enumerate all pairwise non-gap character comparisons at each column and tally substitution pairs.

### Count Substitutions from Alignment

```python
def substitution_counts(alignment):
    # Tally pairwise non-gap residue mismatches across all column-pair comparisons
    ...
```

Full implementation: `examples/substitution_counts.py`. For nucleotide alignments (upper-cased, RNA U mapped to T) it also reports the Ti/Tv ratio over A/C/G/T pairs only, excluding N and other ambiguity codes; for protein alignments it prints no Ti/Tv (Ala/Gly and Cys/Thr are not transitions). The output is raw counts, not a log-odds matrix.

### Built-in Pairwise Substitutions

For pairwise alignments created with `PairwiseAligner`, use the `.substitutions` property:

```python
from Bio.Align import PairwiseAligner
aligner = PairwiseAligner(mode='global', match_score=1, mismatch_score=-1)
seq1, seq2 = (str(r.seq).replace('-', '') for r in alignment[:2])  # ungapped: '-' rows would put '-' in the matrix
print(aligner.align(seq1, seq2)[0].substitutions)
```

### BLOSUM62 Lambda Is Not a Single Number

Different tools calibrate Karlin-Altschul lambda differently (NCBI BLAST tabulates 0.3176; FASTA/SSEARCH recomputes per query; HMMER `phmmer` derives it via Forward calibration; Bio.Align does not implement Karlin-Altschul). Expect ~2% bit-score variation across tools on the same alignment. Always record tool and version when citing bit scores; for borderline (~30 bit) hits this matters.
