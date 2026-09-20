'''Find conserved positions in an alignment'''
# Reference: biopython 1.83+ | Verify API if version differs

from collections import Counter

import numpy as np

from msa_utils import example_path, load_alignment, normalize_alignment


def find_conserved_positions(alignment, threshold=0.8, weights=None):
    '''(column, residue, conservation) for columns whose most common non-gap residue holds >= threshold
    of ALL rows (gap rows count in the denominator, so a column with a gap is never 100% conserved).
    Positions are 0-based. weights: optional per-sequence weights.'''
    alignment = normalize_alignment(alignment)
    weights = np.ones(len(alignment)) if weights is None else np.asarray(weights, dtype=float)
    if len(weights) != len(alignment):
        raise ValueError('weights must have one value per sequence')
    total = weights.sum()
    conserved = []
    for col_idx in range(alignment.get_alignment_length()):
        counts = Counter()
        for char, weight in zip(alignment[:, col_idx], weights):
            counts[char] += weight
        counts.pop('-', None)
        if not counts:
            continue
        most_common_char, most_common_count = counts.most_common(1)[0]
        conservation = most_common_count / total
        if conservation >= threshold - 1e-12:
            conserved.append((col_idx, most_common_char, conservation))
    return conserved


if __name__ == '__main__':
    alignment = load_alignment(example_path('example_alignment.fasta'))

    print('Fully conserved positions (100%):')
    for pos, char, cons in find_conserved_positions(alignment, threshold=1.0):
        print(f'  Position {pos}: {char}')

    print('\nHighly conserved positions (80%+):')
    for pos, char, cons in find_conserved_positions(alignment, threshold=0.8):
        print(f'  Position {pos}: {char} ({cons*100:.0f}%)')
