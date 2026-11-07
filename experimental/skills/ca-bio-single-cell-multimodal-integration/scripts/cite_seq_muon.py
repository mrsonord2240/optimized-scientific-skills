#!/usr/bin/env python3
"""CITE-seq in Python/Muon, CLR-only fallback: QC, per-modality PCA and neighbours, Muon WNN, Leiden, UMAP.

Usage:  python cite_seq_muon.py filtered_feature_bc_matrix.h5 out_prefix [--resolution 0.5]
Input:  10x filtered .h5 with RNA and protein (ADT).
Output: out_prefix.h5mu (obs 'wnn_clusters', obsm X_umap), out_prefix_umap.pdf.
CLR does NOT remove ADT background and Python has no first-party DSB: for denoised protein use
routes/cite-seq-wnn-r.md (DSB) or routes/totalvi.md.
Needs:  muon, scanpy. Checked: scanpy 1.10+, muon (not yet executed after the recut).
"""
import argparse

import muon as mu
import scanpy as sc

ap = argparse.ArgumentParser()
ap.add_argument('h5')
ap.add_argument('out_prefix')
ap.add_argument('--resolution', type=float, default=0.5)
args = ap.parse_args()

mdata = mu.read_10x_h5(args.h5)
rna = mdata.mod['rna']
adt = mdata.mod['prot']

rna.var['mt'] = rna.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(rna, qc_vars=['mt'], percent_top=None, inplace=True)
mu.pp.intersect_obs(mdata)  # keep cells present in both modalities
mask = (rna.obs['n_genes_by_counts'] > 200) & (rna.obs['n_genes_by_counts'] < 5000) & (rna.obs['pct_counts_mt'] < 20)
mdata = mdata[mask].copy()
rna = mdata.mod['rna']  # slicing returns a new MuData: rebind modalities
adt = mdata.mod['prot']

sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
sc.pp.highly_variable_genes(rna, n_top_genes=2000)
sc.pp.scale(rna, max_value=10)
sc.tl.pca(rna, n_comps=30)
sc.pp.neighbors(rna, n_neighbors=15, n_pcs=min(30, rna.obsm['X_pca'].shape[1]))  # Muon WNN combines existing graphs; it does not create them

mu.prot.pp.clr(adt, axis=0)  # the CLR axis convention differs across versions (axis=0 is per feature, Seurat margin=2): check the muon docs
sc.pp.scale(adt, max_value=10)
sc.tl.pca(adt, n_comps=min(18, adt.n_vars - 1))
sc.pp.neighbors(adt, n_neighbors=15, n_pcs=min(18, adt.obsm['X_pca'].shape[1]))

joint_candidates = min(200, mdata.n_obs - 1)  # Muon's default 200 exceeds small experiments
mu.pp.neighbors(mdata, n_neighbors=min(15, joint_candidates), n_bandwidth_neighbors=min(20, joint_candidates),
                n_multineighbors=joint_candidates, key_added='wnn')
sc.tl.leiden(mdata, resolution=args.resolution, key_added='wnn_clusters', neighbors_key='wnn')
mu.tl.umap(mdata, neighbors_key='wnn', random_state=0)  # MuData-aware UMAP; scanpy's AnnData umap cannot read Muon WNN metadata

sc.settings.figdir = './'
sc.settings.autoshow = False
sc.pl.umap(mdata, color='wnn_clusters', save=f'_{args.out_prefix}_umap.pdf')
mdata.write(f'{args.out_prefix}.h5mu')
print(f'wrote {args.out_prefix}.h5mu')
