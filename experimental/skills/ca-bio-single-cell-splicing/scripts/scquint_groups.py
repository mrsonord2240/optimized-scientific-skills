#!/usr/bin/env python3
"""Compare two plate-data cell groups: STARsolo SJ dir, GTF, cell labels, groups and output prefix."""
import argparse
import numpy as np
import pandas as pd
from scquint.data import load_adata_from_starsolo, add_gene_annotation, group_introns
from scquint.differential_splicing import run_differential_splicing

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('sj_dir')
    ap.add_argument('gtf')
    ap.add_argument('labels', help='TSV: first column cell ID; column cell_type')
    ap.add_argument('group1')
    ap.add_argument('group2')
    ap.add_argument('out_prefix')
    a = ap.parse_args()
    adata = load_adata_from_starsolo(a.sj_dir)
    adata = add_gene_annotation(adata, a.gtf)
    if adata.n_vars == 0 or not adata.var['gene_id'].notna().any():
        raise ValueError('zero annotated introns; check chr-prefixed junctions and unprefixed GTF')
    adata = group_introns(adata, by='three_prime')
    labels = pd.read_csv(a.labels, sep='\t', index_col=0)
    if not labels.index.is_unique or 'cell_type' not in labels:
        raise ValueError('labels need unique cell IDs and a cell_type column')
    labels = labels.reindex(adata.obs_names)
    if labels.cell_type.isna().any():
        raise ValueError('missing cell_type for one or more cells')
    adata.obs['cell_type'] = labels.cell_type
    a_idx = np.where(adata.obs.cell_type == a.group1)[0]
    b_idx = np.where(adata.obs.cell_type == a.group2)[0]
    if not len(a_idx) or not len(b_idx):
        raise ValueError('both requested groups must contain cells')
    groups, introns = run_differential_splicing(adata, a_idx, b_idx,
        min_cells_per_intron_group=10, min_total_cells_per_intron=10)
    groups.to_csv(a.out_prefix + '.groups.tsv', sep='\t', index=True)
    introns.to_csv(a.out_prefix + '.introns.tsv', sep='\t', index=True)

if __name__ == '__main__':
    main()
