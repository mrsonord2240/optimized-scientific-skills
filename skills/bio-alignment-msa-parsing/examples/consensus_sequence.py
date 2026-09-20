'''Generate consensus sequence from alignment'''
# Reference: biopython 1.83+ | Verify API if version differs

from collections import Counter

import numpy as np

from msa_utils import example_path, is_nucleotide, load_alignment, normalize_alignment


def consensus_sequence(alignment, threshold=0.5, gap_char='-', ambiguous=None, weights=None):
    '''Most common non-gap character per column when its share of ALL rows (gap rows included)
    is >= threshold, else `ambiguous`. ambiguous=None picks 'N' for nucleotide and 'X' for protein
    ('N' is asparagine in a protein consensus). weights: optional per-sequence weights.'''
    alignment = normalize_alignment(alignment)
    if ambiguous is None:
        ambiguous = 'N' if is_nucleotide(alignment) else 'X'
    weights = np.ones(len(alignment)) if weights is None else np.asarray(weights, dtype=float)
    if len(weights) != len(alignment):
        raise ValueError('weights must have one value per sequence')
    total = weights.sum()
    consensus = []
    for col_idx in range(alignment.get_alignment_length()):
        counts = Counter()
        for char, weight in zip(alignment[:, col_idx], weights):
            counts[char] += weight
        counts.pop('-', None)
        if not counts:
            consensus.append(gap_char)
            continue
        most_common_char, most_common_count = counts.most_common(1)[0]
        consensus.append(most_common_char if most_common_count / total >= threshold else ambiguous)
    return ''.join(consensus)


if __name__ == '__main__':
    alignment = load_alignment(example_path('example_alignment.fasta'))
    print(f'Alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns\n')

    consensus_50 = consensus_sequence(alignment, threshold=0.5)
    print(f'Consensus (50% threshold):\n{consensus_50}\n')

    consensus_70 = consensus_sequence(alignment, threshold=0.7)
    print(f'Consensus (70% threshold):\n{consensus_70}\n')

    consensus_100 = consensus_sequence(alignment, threshold=1.0)
    print(f'Consensus (100% threshold):\n{consensus_100}')
