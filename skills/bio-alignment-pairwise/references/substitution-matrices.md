# Substitution Matrices and Affine Gap Rationale (reference)

## Substitution Matrix Selection

**Goal:** Select the appropriate substitution matrix based on expected sequence divergence.

**Approach:** Match matrix to divergence level. BLOSUM and PAM number in **opposite directions**: higher BLOSUM = closer sequences; higher PAM = more distant sequences.

| Divergence Level | BLOSUM | PAM | When To Use |
|-----------------|--------|-----|-------------|
| Very close (<20% divergence) | BLOSUM80, BLOSUM90 | PAM30 | Recently duplicated genes, strain comparison |
| Moderate | BLOSUM62 (default) | PAM120 | General-purpose, most analyses |
| Distant (>50% divergence) | BLOSUM45, BLOSUM50 | PAM250 | Remote homology detection |

**BLOSUM62 is the universal default** (used by BLAST, most alignment tools). When in doubt, use BLOSUM62. Switch to BLOSUM80 for very similar proteins or BLOSUM45 for distant homologs.

**DNA matrices**: `NUC.4.4` (match=+5, mismatch=-4) accepts IUPAC ambiguity codes but scores partial matches weakly (`R` vs `A` = +1, not +5; check with `print(substitution_matrices.load('NUC.4.4'))` before relying on ambiguity-aware scoring). `HOXD70` is tuned for human-mouse whole-genome alignment from noncoding regions.

```python
from Bio.Align import substitution_matrices
print(substitution_matrices.load())  # List all 30 available matrices

blosum62 = substitution_matrices.load('BLOSUM62')  # General protein (default)
blosum80 = substitution_matrices.load('BLOSUM80')  # Close homologs
blosum45 = substitution_matrices.load('BLOSUM45')  # Distant homologs
nuc44 = substitution_matrices.load('NUC.4.4')      # DNA with IUPAC support
```

### Affine Gap Penalties: Biological Rationale

Gap penalties control how gaps (insertions/deletions) are scored. The **affine model** (`penalty = open + extend * (L-1)`) is almost always preferred over linear because it reflects indel biology: a DNA break introduces the first gap (costly), but extending an existing gap is mechanistically easier (less costly). This models the observation that indels in real sequences tend to occur as single contiguous events.

Typical values with BLOSUM62: `open_gap_score=-12, extend_gap_score=-1` (BLASTP defaults; see SKILL.md Gap Penalties for the convention) or `-11/-1` (EMBOSS). Setting gap open equal to gap extend (linear model) over-penalizes long indels and under-penalizes scattered single-residue gaps, producing biologically unrealistic alignments.
