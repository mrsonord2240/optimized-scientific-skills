# Alignment quality scores (alignment_score and sum of pairs)

Read when scoring an alignment with flat match/mismatch/gap scores or BLOSUM62 sum-of-pairs, or when weighting pairs on an unbalanced dataset. Assumes the normalised `alignment` and the imports from "Required Import and Normalisation" in `SKILL.md`.

## Alignment Quality Metrics

**Goal:** Score alignment quality using sum-of-pairs or simple match/mismatch/gap scoring across all columns.

**Approach:** For each column, score all pairwise residue comparisons and sum across the alignment. Two conventions, both giving a gap/gap pair score 0 (it is not a comparison, and an all-gap column must not change the score): `alignment_score` charges residue/gap pairs a flat `gap` penalty; `sum_of_pairs` skips them. For proteins use the BLOSUM62 form; for DNA, simple match/mismatch.

`examples/alignment_scores.py` holds both functions (importable; `selftest.py` checks them): `python examples/alignment_scores.py [alignment]` prints both scores.

### Overall Alignment Score
```python
from alignment_scores import alignment_score   # examples/alignment_scores.py

score = alignment_score(alignment)
print(f'Alignment score: {score}')
```

### Sum of Pairs Score

Biopython's `substitution_matrices.load('BLOSUM62')` returns a `Bio.Align.substitution_matrices.Array` object (a numpy-backed 2D array indexed by residue characters), not a dict: use `matrix[c1, c2]`. Standard BLOSUM62 includes `B`, `Z`, `X`, and `*`; pairs containing residues outside the matrix alphabet (`U` selenocysteine, `J` Leu/Ile, lower-case letters, `.`) raise `IndexError`; `sum_of_pairs` skips them and warns with the count, so normalise first.

```python
from alignment_scores import sum_of_pairs   # examples/alignment_scores.py

print(f'Sum of pairs (BLOSUM62): {sum_of_pairs(alignment)}')
```

**SP-score is biased on unbalanced datasets.** The implementation gives equal weight to every sequence pair. On phylogenetically structured datasets (e.g. 95 mammals + 5 outgroups), 99% of pairs are mammal-mammal and the SP score reports only mammal-internal alignment quality. MUSCLE and T-Coffee internally compute weighted SP using sequence weights (Henikoff or position-based) so pair contributions are downweighted by cluster redundancy. For SP-as-quality-score on real data, multiply each pair contribution by `weight[i] * weight[j]` from the Henikoff weights in `alignment/msa-parsing` (`examples/henikoff_weights.py`).
