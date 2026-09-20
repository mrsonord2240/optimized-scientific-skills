'''Calculate conservation score at each position'''
# Reference: biopython 1.83+ (checked on 1.88) | Verify API if version differs

from collections import Counter

import numpy as np

from msa_utils import example_path, load_alignment


def column_conservation(alignment, col_idx, ignore_gaps=True, min_occupancy=0.5):
    '''Frequency of the most common residue; NaN for empty columns and (with ignore_gaps) for columns
    where fewer than min_occupancy of the sequences have a residue (2 residues of 73 are not 100% conserved).'''
    full = alignment[:, col_idx]
    column = full.replace('-', '') if ignore_gaps else full
    if not column or (ignore_gaps and len(column) / len(full) < min_occupancy):
        return float('nan')
    return Counter(column).most_common(1)[0][1] / len(column)


def conservation_profile(alignment, window=10, min_occupancy=0.5):
    '''Centred moving average over columns i-window//2 .. i+window//2 (inclusive), NaN columns skipped.'''
    length = alignment.get_alignment_length()
    scores = np.array([column_conservation(alignment, i, min_occupancy=min_occupancy) for i in range(length)])
    half = window // 2
    profile = []
    for i in range(length):
        chunk = scores[max(0, i - half):i + half + 1]
        chunk = chunk[~np.isnan(chunk)]
        profile.append(chunk.mean() if chunk.size else float('nan'))
    return profile


if __name__ == '__main__':
    alignment = load_alignment(example_path('example_protein.fasta'))
    print(f'Alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns\n')

    conservation_scores = [column_conservation(alignment, i) for i in range(alignment.get_alignment_length())]

    print('Conservation profile (NaN = fewer than 50% of sequences have a residue):')
    for i, score in enumerate(conservation_scores):
        bar = '' if np.isnan(score) else '#' * int(score * 20)
        print(f'  {i:3d}: {score*100:5.1f}% {bar}')

    used = [s for s in conservation_scores if not np.isnan(s)]
    print(f'\nAverage conservation: {np.mean(used)*100:.1f}% over {len(used)} of {len(conservation_scores)} columns')

    fully_conserved = sum(1 for s in used if s == 1.0)
    print(f'Fully conserved columns: {fully_conserved} ({fully_conserved/len(conservation_scores)*100:.1f}%)')
