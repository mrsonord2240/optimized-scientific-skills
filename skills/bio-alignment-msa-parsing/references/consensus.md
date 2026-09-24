# Consensus Sequence

Read when the request is for a consensus sequence. Runs `examples/consensus_sequence.py`; standalone: `python examples/consensus_sequence.py [alignment_file]`.

## Consensus Sequence

**"Get consensus sequence"** -> Derive a single representative sequence from an MSA based on majority-rule voting at each column.

**Goal:** Generate a consensus sequence from the alignment using a frequency threshold.

**Approach:** At each column, select the most common non-gap character if its share of all rows (gap rows included) reaches the threshold; otherwise mark as ambiguous. The placeholder is alphabet-aware: `N` for nucleotide, `X` for protein (`N` is asparagine and cannot be told apart from a real Asn column).

### Simple Majority Consensus
```python
from consensus_sequence import consensus_sequence  # examples/consensus_sequence.py (sys.path as in SKILL.md)

consensus = consensus_sequence(alignment, threshold=0.5)
```

Signature: `consensus_sequence(alignment, threshold=0.5, gap_char='-', ambiguous=None, weights=None)`. `ambiguous=None` picks `N` for nucleotide and `X` for protein (alphabet test `is_nucleotide()` in `examples/msa_utils.py`: at least 90% of the non-gap characters are ACGTUN); an all-gap column gives `gap_char`; `weights` of the wrong length or summing to zero raise `ValueError`.

### Note on Bio.Align.AlignInfo
`AlignInfo.SummaryInfo` keeps only `get_column` in Biopython 1.88 (`dumb_consensus`, `gap_consensus`, `pos_specific_score_matrix`, `information_content` were removed and raise `AttributeError`). Use `consensus_sequence()` from `examples/consensus_sequence.py`.
