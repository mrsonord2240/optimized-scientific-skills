# Directed fate mapping with CellRank 2: PseudotimeKernel + ConnectivityKernel -> GPCCA.
# Input : .h5ad with a pseudotime column in .obs (default dpt_pseudotime) and a cluster column
#         (default leiden), neighbors already computed.
# Output: prints macrostates, terminal/initial states and per-cluster mean fate-probability entropy;
#         with --out, writes the AnnData carrying fate probabilities.
# Usage : python cellrank_fate.py adata.h5ad [--time-key dpt_pseudotime] [--cluster-key leiden]
#                                 [--n-states 10] [--w-pseudotime 0.8] [--out fate.h5ad]
# Checked on cellrank 2.3.3. The __main__ guard is required on Windows (multiprocessing progress bar).
import argparse

import numpy as np
import pandas as pd
import scanpy as sc
import cellrank as cr


def run_fate_mapping(adata, time_key='dpt_pseudotime', cluster_key='leiden', n_states=10, w_pseudotime=0.8):
    pk = cr.kernels.PseudotimeKernel(adata, time_key=time_key).compute_transition_matrix(n_jobs=1)
    ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
    combined = w_pseudotime * pk + (1 - w_pseudotime) * ck   # weights are a researcher choice; sweep them

    g = cr.estimators.GPCCA(combined)
    g.compute_macrostates(n_states=n_states, cluster_key=cluster_key)   # n_states from the Schur/eigenvalue spectral gap
    g.predict_terminal_states(method='stability')
    g.predict_initial_states(n_states=1, allow_overlap=True)   # without allow_overlap, real branching data can raise
                                                                 # ValueError: N cells overlapped between initial/terminal states
    g.compute_fate_probabilities(n_jobs=1)
    g.compute_lineage_drivers()
    return g


def fate_entropy(g, adata):
    fp = g.fate_probabilities
    df = pd.DataFrame(np.asarray(fp.X), index=adata.obs_names, columns=fp.names).clip(lower=1e-12)
    return -(df * np.log(df)).sum(axis=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('h5ad')
    ap.add_argument('--time-key', default='dpt_pseudotime')
    ap.add_argument('--cluster-key', default='leiden')
    ap.add_argument('--n-states', type=int, default=10)
    ap.add_argument('--w-pseudotime', type=float, default=0.8)
    ap.add_argument('--out')
    a = ap.parse_args()

    adata = sc.read_h5ad(a.h5ad)
    g = run_fate_mapping(adata, a.time_key, a.cluster_key, a.n_states, a.w_pseudotime)
    print('macrostates:', list(g.macrostates.cat.categories))
    print('terminal_states:', list(g.terminal_states.cat.categories.dropna()))
    print('initial_states:', list(g.initial_states.cat.categories.dropna()))
    print('mean fate-probability entropy per cluster:')
    print(fate_entropy(g, adata).groupby(adata.obs[a.cluster_key].values).mean().round(3).to_string())
    if a.out:
        adata.write(a.out)


if __name__ == '__main__':
    main()
