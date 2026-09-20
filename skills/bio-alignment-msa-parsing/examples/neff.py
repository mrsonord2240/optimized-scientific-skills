'''Compute effective sequence number (Neff) by similarity-based clustering.

For each sequence count the sequences (itself included) sharing >= identity_threshold of the
positions where both have a residue; each sequence gets weight 1 / count and Neff is the sum of
the weights. Identity 0.62 is the hmmbuild --wblosum default (BLOSUM62's clustering threshold);
0.80 is a common DCA convention for nucleotide/deep protein MSAs.

This is an Neff (a count of non-redundant sequences). pyhmmer compute_weights() returns per-sequence
weights that sum to N instead, and hmmbuild's eff_nseq is an entropy-weighted quantity; the three
differ by up to 10x (see SKILL.md).
'''
# Reference: biopython 1.83+, numpy 1.26+ | Verify API if version differs

import numpy as np

from msa_utils import example_path, load_alignment, to_array

# Single rule used across this Skill for MI-APC / DCA: Neff/L must exceed this value.
MIN_NEFF_PER_L = 1.0


def neff(alignment, identity_threshold=0.62):
    seq_array = to_array(alignment)
    n = len(seq_array)
    cluster_size = np.ones(n)
    for i in range(n):
        for j in range(i + 1, n):
            mask = (seq_array[i] != '-') & (seq_array[j] != '-')
            length = mask.sum()
            if not length:
                continue
            shared = ((seq_array[i] == seq_array[j]) & mask).sum()
            if shared / length >= identity_threshold:
                cluster_size[i] += 1
                cluster_size[j] += 1
    weights = 1.0 / cluster_size
    return weights.sum()


if __name__ == '__main__':
    alignment = load_alignment(example_path('example_alignment.fasta'))
    n_eff_protein = neff(alignment, identity_threshold=0.62)
    n_eff_nucleotide = neff(alignment, identity_threshold=0.80)
    length = alignment.get_alignment_length()

    print(f'Sequences: {len(alignment)}')
    print(f'Length: {length}')
    print(f'Neff (62% threshold, protein convention): {n_eff_protein:.2f}')
    print(f'Neff/L: {n_eff_protein / length:.3f}')
    print(f'Neff (80% threshold, nucleotide convention): {n_eff_nucleotide:.2f}')
    verdict = 'meets' if n_eff_protein / length > MIN_NEFF_PER_L else 'is below'
    print(f'Neff/L {verdict} the Skill rule of thumb for MI-APC / DCA (Neff/L > {MIN_NEFF_PER_L:g}).')
