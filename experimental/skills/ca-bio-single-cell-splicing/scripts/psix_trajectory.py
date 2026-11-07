#!/usr/bin/env python3
"""Run the source Psix trajectory workflow: cells.h5ad psi.tsv mrna.tsv results.tsv."""
import argparse
import pandas as pd
import scanpy as sc
import psix

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('cells_h5ad')
    ap.add_argument('psi_table')
    ap.add_argument('mrna_table')
    ap.add_argument('output')
    ap.add_argument('--jobs', type=int, default=4)
    a = ap.parse_args()
    adata = sc.read_h5ad(a.cells_h5ad)
    if 'X_pca' not in adata.obsm:
        raise ValueError('cells.h5ad must contain X_pca')
    latent = pd.DataFrame(adata.obsm['X_pca'][:, :20], index=adata.obs_names)
    model = psix.Psix(psi_table=a.psi_table, mrna_table=a.mrna_table)
    model.run_psix(latent=latent, n_jobs=a.jobs)
    model.psix_results.to_csv(a.output, sep='\t', index=True)

if __name__ == '__main__':
    main()
