'''Calculate Shannon entropy and Kullback-Leibler information content per column.'''
# Reference: biopython 1.83+ (checked on 1.88) | Verify API if version differs

import math
from collections import Counter

from msa_utils import check_alphabet, example_path, load_alignment, pick_background


def shannon_entropy(column, ignore_gaps=True):
    if ignore_gaps:
        column = column.replace('-', '')
    if not column:
        return 0.0
    counts = Counter(column)
    total = len(column)
    entropy = 0.0
    for count in counts.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy


def information_content(column, background, ignore_gaps=True):
    '''KL-divergence IC against a background distribution, in bits.

    Use ROBINSON_BACKGROUND for protein or DNA_UNIFORM for nucleotide alignments. Letters outside the
    background (X, B, Z, U, N, ...) are dropped and the column renormalised over the rest; an empty
    column returns 0.0. Input must be normalised (upper case, '-' gaps): run normalize_alignment() first.
    '''
    letters = [r for r in column if r in background]
    if not letters:
        return 0.0
    total = len(letters)
    return sum((c / total) * math.log2((c / total) / background[r]) for r, c in Counter(letters).items())


if __name__ == '__main__':
    alignment = load_alignment(example_path('example_protein.fasta'))
    print(f'Alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns\n')

    background, alphabet_label = pick_background(alignment)
    print(f'Treating as {alphabet_label}')
    check_alphabet(alignment, background, 'information content')
    print()

    print('Position   Entropy (bits)   IC (bits)')
    print('-' * 45)
    entropies = []
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        ent = shannon_entropy(column)
        ic = information_content(column, background)
        entropies.append(ent)
        print(f'{i:5d}      {ent:8.3f}         {ic:8.3f}')

    avg_entropy = sum(entropies) / len(entropies)
    print(f'\nAverage entropy: {avg_entropy:.3f} bits')
    print(f'Maximum possible entropy: {math.log2(len(background)):.3f} bits')
