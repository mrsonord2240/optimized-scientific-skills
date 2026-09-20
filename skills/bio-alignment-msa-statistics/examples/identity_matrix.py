'''Pairwise identity matrix for an alignment (PID1-PID4, Raghava & Barton 2006).

Reference: BioPython 1.83+ (checked on 1.88), numpy 1.26+ | Verify API if version differs

Usage: python identity_matrix.py [alignment.fasta] [pid1|pid2|pid3|pid4]   (default: shipped example, pid4)

PID1 = identical / (aligned pairs + internal gap columns): only the span from the pair's first to its last
aligned (both-residue) column counts, so terminal overhangs and unaligned flanks are excluded
(pwalign::pid type PID1).
PID2 = identical / aligned pairs. PID3 = identical / shorter ungapped length.
PID4 = identical / mean ungapped length (default; recommended for evolutionary analyses).
Undefined identities (an all-gap sequence, or no aligned pair for PID1/PID2) are NaN, never 0.

Vectorized via numpy broadcasting; O(N^2 * L) but with numpy-internal loops.
Above ~10k sequences, switch to k-mer-based distance estimation (e.g. mash).
'''
import sys

import numpy as np

from msa_utils import example_path, load_alignment

METHODS = ('pid1', 'pid2', 'pid3', 'pid4')


def pairwise_identity(seq1, seq2, method='pid4'):
    '''Identity of two normalised aligned rows ('-' gaps, upper case); NaN when undefined.'''
    if not seq1.replace('-', '') or not seq2.replace('-', ''):
        return float('nan')  # an all-gap sequence has no identity to anything
    both = [i for i, (a, b) in enumerate(zip(seq1, seq2)) if a != '-' and b != '-']
    matches = sum(seq1[i] == seq2[i] for i in both)
    if method == 'pid1':
        span = range(both[0], both[-1] + 1) if both else range(0)
        denom = sum(seq1[i] != '-' or seq2[i] != '-' for i in span)
    elif method == 'pid2':
        denom = len(both)
    elif method == 'pid3':
        denom = min(len(seq1.replace('-', '')), len(seq2.replace('-', '')))
    elif method == 'pid4':
        denom = (len(seq1.replace('-', '')) + len(seq2.replace('-', ''))) / 2
    else:
        raise ValueError(f'method must be one of {METHODS}, got {method!r}')
    return matches / denom if denom > 0 else float('nan')


def identity_matrix_vectorized(alignment, method='pid4'):
    '''N x N identity matrix (fractions) of a normalised alignment; NaN where undefined.'''
    if method not in METHODS:
        raise ValueError(f'method must be one of {METHODS}, got {method!r}')
    arr = np.array([list(str(r.seq)) for r in alignment])
    residue = arr != '-'
    n, length = residue.shape
    columns = np.arange(length)
    lengths = residue.sum(axis=1)
    matrix = np.full((n, n), np.nan)
    for i in range(n):
        both = residue[i] & residue
        matches = ((arr[i] == arr) & both).sum(axis=1)
        if method == 'pid1':
            first = both.argmax(axis=1)
            last = length - 1 - both[:, ::-1].argmax(axis=1)
            in_span = (columns >= first[:, None]) & (columns <= last[:, None])
            denom = ((residue[i] | residue) & in_span).sum(axis=1) * both.any(axis=1)
        elif method == 'pid2':
            denom = both.sum(axis=1)
        elif method == 'pid3':
            denom = np.minimum(lengths[i], lengths)
        else:
            denom = (lengths[i] + lengths) / 2
        denom = np.where((lengths[i] > 0) & (lengths > 0), denom, 0)
        matrix[i] = np.where(denom > 0, matches / np.maximum(denom, 1), np.nan)
    return matrix


def average_identity(matrix):
    '''Mean off-diagonal identity ignoring NaN, and how many pairs were undefined.'''
    off = matrix[~np.eye(len(matrix), dtype=bool)]
    finite = off[np.isfinite(off)]
    return (finite.mean() if finite.size else float('nan')), int(off.size - finite.size)


if __name__ == '__main__':
    method = sys.argv[2] if len(sys.argv) > 2 else 'pid4'
    alignment = load_alignment(example_path('example_protein.fasta'))
    seq_ids = [r.id for r in alignment]

    matrix = identity_matrix_vectorized(alignment, method)

    print(f'Pairwise Identity Matrix (%, {method.upper()}):')
    print(f'{"":>12}', ' '.join(f'{s[:8]:>8}' for s in seq_ids))
    for i, row in enumerate(matrix):
        print(f'{seq_ids[i][:12]:>12}', ' '.join(f'{v*100:>7.1f}%' for v in row))

    avg_identity, undefined = average_identity(matrix)
    print(f'\nAverage pairwise identity: {avg_identity*100:.1f}%  ({undefined} undefined pairs ignored)')
