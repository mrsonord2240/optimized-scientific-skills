'''Target-decoy FDR and q-values from a concatenated-search PSM table.

Demonstrates the load-bearing idea: a raw PSM score is meaningless in isolation;
the actionable number is a list-level q-value (or a per-PSM PEP). This uses the
CONCATENATED target-decoy competition estimator FDR = (decoys + 1)/targets on
one best hit per spectrum. Separate target/decoy searches would instead need
pi0 * decoys/targets (Kall 2008) or the mix-max estimator (Keich 2015).
Self-contained: builds a synthetic table so no input files are needed.'''
# Reference: numpy 1.26+, pandas 2.2+ | Verify API if version differs
import numpy as np
import pandas as pd

DECOY_PREFIXES = ('decoy_', 'rev_', 'xxx_')   # compared lower-cased; Sage and FragPipe write rev_
TARGET_FDR = 0.01   # 1% list-level FDR, the community standard for peptide IDs


def build_demo_table(n_true=1000, n_null=3000, seed=0):
    '''Concatenated search: every spectrum gets one best target and one best decoy
    candidate and keeps the higher (target-decoy competition). True spectra have a
    high-scoring correct target; null spectra draw target and decoy from one null.'''
    rng = np.random.default_rng(seed)
    n = n_true + n_null
    target = np.concatenate([rng.normal(3.5, 0.8, n_true), rng.normal(1.0, 0.8, n_null)])
    decoy = rng.normal(1.0, 0.8, n)
    decoy_wins = decoy > target
    proteins = np.where(decoy_wins, f'{DECOY_PREFIXES[0]}prot', 'TARGET')
    return pd.DataFrame({'scan': np.arange(n), 'protein': proteins,
                         'score': np.maximum(target, decoy),
                         'is_true': (~decoy_wins) & (np.arange(n) < n_true)})


def add_qvalues(psms):
    '''Concatenated competition: best hit per spectrum, rank by score,
    FDR = (cumulative decoys + 1) / targets, then take the running minimum from
    the bottom to make q-values monotone.'''
    psms = psms.copy()
    psms['is_decoy'] = psms['protein'].str.lower().str.startswith(DECOY_PREFIXES)
    if not psms['is_decoy'].any():
        raise ValueError('no decoy PSMs recognised: check the decoy prefix')
    psms = psms.sort_values('score', ascending=False).drop_duplicates('scan').reset_index(drop=True)
    targets = (~psms['is_decoy']).cumsum()
    decoys = psms['is_decoy'].cumsum()
    psms['fdr'] = (decoys + 1) / targets.clip(lower=1)
    psms['qvalue'] = psms['fdr'][::-1].cummin()[::-1]
    return psms


def add_pep(psms, bandwidth=0.3):
    '''PEP (local FDR) is per-PSM: the decoy density over total density at a score.
    Estimated here by a simple kernel-smoothed decoy fraction in a score window.'''
    psms = psms.copy()
    s = psms['score'].to_numpy()
    is_decoy = psms['is_decoy'].to_numpy().astype(float)
    pep = np.empty(len(s))
    for i, score in enumerate(s):
        w = np.exp(-((s - score) ** 2) / (2 * bandwidth ** 2))
        decoy_density = (is_decoy * w).sum()
        target_density = ((1 - is_decoy) * w).sum()
        pep[i] = min(1.0, decoy_density / max(target_density, 1e-9))   # local false-target rate, concatenated competition (no x2)
    psms['pep'] = pep
    return psms


def main():
    psms = add_qvalues(build_demo_table())
    psms = add_pep(psms)
    n_t = (~psms['is_decoy']).sum()
    n_d = psms['is_decoy'].sum()
    print(f'Targets: {n_t}, Decoys: {n_d} (concatenated 1:1 search)')

    kept = psms[(psms['qvalue'] <= TARGET_FDR) & (~psms['is_decoy'])]
    print(f'Target PSMs at q <= {TARGET_FDR}: {len(kept)}')
    print(f'Worst PEP inside the {TARGET_FDR:.0%}-FDR list: {kept["pep"].max():.3f}')

    strict = psms[(psms['pep'] <= TARGET_FDR) & (~psms['is_decoy'])]
    print(f'Target PSMs at PEP <= {TARGET_FDR}: {len(strict)} '
          f'(per-PSM cutoff is far stricter than the same q-value)')


if __name__ == '__main__':
    main()
