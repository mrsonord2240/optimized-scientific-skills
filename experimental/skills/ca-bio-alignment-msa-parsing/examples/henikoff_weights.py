'''Compute Henikoff position-based sequence weights for an MSA.

Reference: Henikoff S, Henikoff JG. 1994. Position-based sequence weights. JMB 243:574-578.
Each column contributes 1 / (k * n) per sequence, where k is the number of distinct residues
and n is the count of the sequence's residue at that column.

Gap handling: columns containing any gap are skipped entirely (weights sum to 1). This is NOT
what HMMER/Easel does: pyhmmer's compute_weights('pb') ignores gaps column by column, uses only
"consensus" columns (>= 50% residues by default), divides by each sequence's residue count and
rescales to sum to N. On the Pfam PF00042 seed the two agree in rank (Spearman 0.909) but not in
value (max absolute difference 0.0091 after rescaling to sum 1). Sequences whose residues sit only in
gappy columns get weight 0 here; for gappy or fragmentary alignments use pyhmmer (see SKILL.md).
'''
# Reference: biopython 1.83+, numpy 1.26+ | Verify API if version differs

import numpy as np

from msa_utils import example_path, load_alignment, to_array


def henikoff_weights(alignment):
    seq_array = to_array(alignment)
    weights = np.zeros(len(seq_array))
    used_columns = 0
    for col_idx in range(seq_array.shape[1]):
        residues, inverse, counts = np.unique(seq_array[:, col_idx], return_inverse=True, return_counts=True)
        if '-' in residues:
            continue
        used_columns += 1
        weights += 1.0 / (len(residues) * counts[inverse])
    if not used_columns:
        raise ValueError('every column contains a gap, so no column can be weighted; '
                         'use pyhmmer compute_weights(method="pb") or trim gappy columns first')
    return weights / weights.sum()


if __name__ == '__main__':
    alignment = load_alignment(example_path('example_alignment.fasta'))
    weights = henikoff_weights(alignment)
    for record, weight in zip(alignment, weights):
        print(f'{record.id}: {weight:.4f}')
    print(f'\nTotal weight: {weights.sum():.4f}  Kish effective sample size (1/sum w^2, not Neff): {1 / (weights ** 2).sum():.2f}')
