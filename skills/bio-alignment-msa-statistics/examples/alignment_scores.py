'''Alignment quality scores from SKILL.md: flat-gap alignment_score and BLOSUM62 sum_of_pairs.

Both give a gap/gap pair the score 0 (not a comparison, so an all-gap column does not change the score).
alignment_score charges residue/gap pairs a flat `gap` penalty; sum_of_pairs skips them and counts the
residue pairs outside the matrix alphabet (U, J, lower case, '.'), so normalise first.
'''
# Reference: biopython 1.83+ (checked on 1.88) | Verify API if version differs

import sys

from Bio.Align import substitution_matrices

from msa_utils import example_path, load_alignment


def alignment_score(alignment, match=1, mismatch=-1, gap=-2):
    total_score = 0
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        for i, c1 in enumerate(column):
            for c2 in column[i+1:]:
                if c1 == '-' and c2 == '-':
                    continue                      # gap/gap is not a comparison
                elif c1 == '-' or c2 == '-':
                    total_score += gap
                elif c1 == c2:
                    total_score += match
                else:
                    total_score += mismatch
    return total_score


def sum_of_pairs(alignment, substitution_matrix=None):
    if substitution_matrix is None:
        substitution_matrix = substitution_matrices.load('BLOSUM62')

    total, skipped = 0.0, 0
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        for i, c1 in enumerate(column):
            for c2 in column[i+1:]:
                if c1 == '-' or c2 == '-':
                    continue
                try:
                    total += substitution_matrix[c1, c2]
                except (KeyError, IndexError):
                    skipped += 1
    if skipped:
        print(f'WARNING: {skipped} residue pairs outside the matrix alphabet were skipped', file=sys.stderr)
    return total


if __name__ == '__main__':
    alignment = load_alignment(example_path('example_protein.fasta'))
    print(f'Alignment score (match 1, mismatch -1, gap -2): {alignment_score(alignment)}')
    print(f'Sum of pairs (BLOSUM62): {sum_of_pairs(alignment)}')
