#!/usr/bin/env python3
"""MOFA+ shared and modality-specific factors on a MuData via muon.

Usage:  python mofa_muon.py in.h5mu out_prefix [--n-factors 15]
Input:  .h5mu with one AnnData per modality, each already in the representation its likelihood needs
        (scaled RNA: gaussian; binarized ATAC: bernoulli; counts: poisson). Likelihoods must match the data.
Output: out_prefix.h5mu with obsm['X_mofa'], and out_prefix.hdf5 (the MOFA model, loadings, variance explained).
MOFA+ is for interpretation (shared vs specific axes of variation), not clustering or denoising, and does not
correct batch within a view.
Needs:  muon, mudata, mofapy2. Checked: not yet executed after the recut.
"""
import argparse

import mudata as md
import muon as mu

ap = argparse.ArgumentParser()
ap.add_argument('h5mu')
ap.add_argument('out_prefix')
ap.add_argument('--n-factors', type=int, default=15)
args = ap.parse_args()

mdata = md.read_h5mu(args.h5mu)
mu.tl.mofa(mdata, n_factors=args.n_factors, outfile=f'{args.out_prefix}.hdf5')  # writes mdata.obsm['X_mofa']
mdata.write_h5mu(f'{args.out_prefix}.h5mu')
print(f"X_mofa {mdata.obsm['X_mofa'].shape}")
