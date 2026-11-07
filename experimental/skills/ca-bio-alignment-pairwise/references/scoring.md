# Scoring

To inspect multiple optimal alignments, use `itertools.islice(aligner.align(target, query), limit)`. Avoid `len()` on the lazy result, which can overflow when many alignments tie. The command-line wrapper emits the first optimum.



### Gap Penalties: Set Them Explicitly, and Mind the Convention

`PairwiseAligner()` with no arguments (Biopython 1.88) uses match_score=1, mismatch_score=0 and open_gap_score=extend_gap_score=-1 (`print(aligner)` shows the resolved set). These defaults are tuned for nothing: combined with BLOSUM62, gaps are nearly free and alignments fragment (HBA vs HBB: 55 gap positions in 37 aligned blocks with the defaults, versus 9 gap positions in 5 blocks at open -11 / extend -1). **Always specify gap penalties explicitly.**

Tools disagree on what "open" means, so the same numbers give different scores:

| Convention | Cost of a gap of length k | Used by |
|---|---|---|
| open + (k-1) * extend | open covers the first gap position | Biopython `open_gap_score`, EMBOSS `needle`/`water`, parasail |
| open + k * extend | open is charged on top of the extension | BLAST+ `-gapopen`, R pwalign `gapOpening` |

Conversion: Biopython `open_gap_score = -(BLAST gapopen + gapextend)`, `extend_gap_score = -gapextend`. So **BLASTP defaults (BLOSUM62, 11/1) are `open_gap_score=-12, extend_gap_score=-1`**, and `-11/-1` is EMBOSS 11/1 (pwalign `gapOpening=10, gapExtension=1`). EMBOSS defaults 10/0.5 are `-10/-0.5`.

Verified on HBA_HUMAN vs HBB_HUMAN, BLOSUM62, local: Biopython -12/-1 = 285 = `blastp -comp_based_stats 0` raw score (same HSP, query 3-141) = pwalign `gapOpening=11, gapExtension=1`; Biopython -11/-1 = 288 = EMBOSS `water -gapopen 11 -gapextend 1` = pwalign `gapOpening=10, gapExtension=1`. Global: -11/-1 = 286 = `needle` 11/1. BLASTP's own default composition-based statistics change the reported score (286 here); use `-comp_based_stats 0` to compare raw scores.

EMBOSS command lines for the same numbers (6.6.0, one sequence per FASTA file; `-auto` suppresses prompts; the matrix file is `EBLOSUM62`, not `BLOSUM62`; DNA uses the default `EDNAFULL`, gaps 10/0.5):

```bash
needle -auto -asequence a.fa -bsequence b.fa -gapopen 11 -gapextend 1 -datafile EBLOSUM62 -outfile needle.txt   # global; "# Score:" line
water  -auto -asequence a.fa -bsequence b.fa -gapopen 11 -gapextend 1 -datafile EBLOSUM62 -outfile water.txt    # local
```

### DNA/RNA Alignment
```python
aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
```

### Protein Alignment
```python
from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load('BLOSUM62')
aligner = PairwiseAligner(mode='global', substitution_matrix=blosum62, open_gap_score=-11, extend_gap_score=-1)  # EMBOSS 11/1; BLASTP-equivalent is -12/-1; other matrices: references/substitution-matrices.md
```

### Local Alignment (Find Best Region)
```python
aligner = PairwiseAligner(mode='local', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
```

### Semiglobal (Overlap/Fragment Alignment)
```python
# Free end gaps on BOTH sequences (order-independent) -- overlap detection between two reads,
# or a fragment inside a reference
aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
aligner.end_gap_score = 0.0
alignment = aligner.align(reference, fragment)[0]   # 20-nt fragment in a 620-nt reference: score 40, target span [[300, 320]]

# Free end gaps on the query (second argument) only. Argument ORDER matters: the fragment must be
# the second argument (align(fragment, reference) scored -279 instead of 40).
aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
aligner.open_left_deletion_score = aligner.extend_left_deletion_score = 0
aligner.open_right_deletion_score = aligner.extend_right_deletion_score = 0
alignment = aligner.align(reference, fragment)[0]
```
Older Biopython named these `query_left_open_gap_score`, `query_left_extend_gap_score`, `query_right_open_gap_score`, `query_right_extend_gap_score` (still accepted on 1.88, with a DeprecationWarning).
