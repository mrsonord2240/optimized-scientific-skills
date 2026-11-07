'''Mutual information with average-product correction (MI-APC) for coevolution detection.

Dunn SD, Wahl LM, Gloor GB. 2008. Mutual information without the influence of phylogeny
or entropy dramatically improves residue contact prediction. Bioinformatics 24:333-340.
APC removes background signal from per-column entropy and phylogenetic structure.
The APC term excludes the i==j diagonal so the per-column mean reflects only
off-diagonal MI contributions, matching the Dunn et al definition.

Guard (the Skill's own rule): APC is applied only when L > 100 and Neff/L > 1. Otherwise
mi_matrix_apc() warns and returns raw MI, and the ranking is only meaningful if it beats the
column-shuffled null printed by main. On the Pfam PF00042 seed (L=141, Neff/L=0.47) the top
score sits below the shuffled null and none of the top pairs are contacts in PDB 1MBN.
'''
# Reference: biopython 1.83+, numpy 1.26+ | Verify API if version differs

import warnings

import numpy as np

from msa_utils import example_path, load_alignment, to_array
from neff import MIN_NEFF_PER_L, neff

APC_MIN_LENGTH = 100


def joint_counts(col_i, col_j):
    pairs = np.array([col_i, col_j])
    unique_i = np.unique(col_i)
    unique_j = np.unique(col_j)
    table = np.zeros((len(unique_i), len(unique_j)))
    idx_i = {r: k for k, r in enumerate(unique_i)}
    idx_j = {r: k for k, r in enumerate(unique_j)}
    for a, b in pairs.T:
        table[idx_i[a], idx_j[b]] += 1
    return table


def mi_from_array(seq_array, min_pairs=20):
    '''Raw pairwise MI (bits) between columns of a (n_seqs, n_cols) character array; gap rows are
    dropped per pair and pairs with fewer than min_pairs shared rows are left at 0.
    Returns (mi, valid) where valid marks the pairs that were scored.'''
    n_cols = seq_array.shape[1]
    mi = np.zeros((n_cols, n_cols))
    valid = np.zeros((n_cols, n_cols), dtype=bool)
    for i in range(n_cols):
        for j in range(i + 1, n_cols):
            col_i = seq_array[:, i]
            col_j = seq_array[:, j]
            mask = (col_i != '-') & (col_j != '-')
            if mask.sum() < min_pairs:
                continue
            table = joint_counts(col_i[mask], col_j[mask])
            joint_p = table / table.sum()
            margin_i = joint_p.sum(axis=1, keepdims=True)
            margin_j = joint_p.sum(axis=0, keepdims=True)
            with np.errstate(divide='ignore', invalid='ignore'):
                pmi = np.where(joint_p > 0, joint_p * np.log2(joint_p / (margin_i * margin_j)), 0)
            mi[i, j] = mi[j, i] = pmi.sum()
            valid[i, j] = valid[j, i] = True
    return mi, valid


def apc_correction(mi, valid):
    n_cols = mi.shape[0]
    off_diag = ~np.eye(n_cols, dtype=bool) & valid
    column_means = np.zeros(n_cols)
    for k in range(n_cols):
        row_mask = off_diag[k]
        if row_mask.any():
            column_means[k] = mi[k, row_mask].mean()
    overall_mean = mi[off_diag].mean() if off_diag.any() else 0.0
    return np.outer(column_means, column_means) / overall_mean if overall_mean > 0 else np.zeros_like(mi)


def apc_applicable(alignment):
    '''(ok, neff_per_l, message): the Skill's guard, L > 100 and Neff/L > 1.'''
    length = alignment.get_alignment_length()
    neff_per_l = neff(alignment) / length
    ok = length > APC_MIN_LENGTH and neff_per_l > MIN_NEFF_PER_L
    message = (f'L={length}, Neff/L={neff_per_l:.2f}: APC needs L > {APC_MIN_LENGTH} and '
               f'Neff/L > {MIN_NEFF_PER_L:g}')
    return ok, neff_per_l, message


def mi_matrix_apc(alignment, min_pairs=20, force=False):
    '''MI - APC when the guard passes (or force=True); otherwise warn and return raw MI.'''
    ok, _, message = apc_applicable(alignment)
    mi, valid = mi_from_array(to_array(alignment), min_pairs)
    if ok or force:
        return mi - apc_correction(mi, valid)
    warnings.warn(f'{message}; returning RAW MI, not MI-APC', stacklevel=2)
    return mi


def column_shuffled_null(alignment, apc=True, min_pairs=20, seed=0, n_shuffles=5):
    '''Best pair score over n_shuffles alignments whose columns were each shuffled independently
    (keeps column entropy, destroys covariation), scored the same way as the real alignment
    (apc=True for MI-APC, False for raw MI). Real pairs must beat this.'''
    rng = np.random.default_rng(seed)
    best = -np.inf
    for _ in range(n_shuffles):
        shuffled = np.array(to_array(alignment))
        for k in range(shuffled.shape[1]):
            shuffled[:, k] = rng.permutation(shuffled[:, k])
        mi, valid = mi_from_array(shuffled, min_pairs)
        score = mi - apc_correction(mi, valid) if apc else mi
        best = max(best, score[np.triu_indices(score.shape[0], k=1)].max())
    return best


if __name__ == '__main__':
    alignment = load_alignment(example_path('example_alignment.fasta'))
    ok, neff_per_l, message = apc_applicable(alignment)
    label = 'MI-APC' if ok else 'raw MI'
    min_pairs = min(20, len(alignment))  # pairs need this many shared non-gap rows
    if not ok:
        print(f'WARNING: {message}. Ranking raw MI; treat it as noise unless it beats the shuffled null.')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        scores = mi_matrix_apc(alignment, min_pairs=min_pairs)
    null_max = column_shuffled_null(alignment, apc=ok, min_pairs=min_pairs)

    upper = np.triu_indices(scores.shape[0], k=1)
    top_pairs = sorted(zip(scores[upper], upper[0], upper[1]), reverse=True)[:20]
    print(f'Top 20 column pairs ({label}, bits); column-shuffled null max = {null_max:.3f}:')
    for score, i, j in top_pairs:
        flag = '' if score > null_max else '  (<= null)'
        print(f'  {i:4d}-{j:4d}: {score:6.3f}{flag}')
    n_above = sum(1 for s, _, _ in top_pairs if s > null_max)
    print(f'{n_above} of the top 20 exceed the null.')
