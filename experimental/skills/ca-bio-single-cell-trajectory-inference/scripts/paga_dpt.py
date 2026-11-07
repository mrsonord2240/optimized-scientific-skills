"""Compute PAGA and rooted DPT on an existing clustered AnnData."""
import argparse
import numpy as np
import scanpy as sc

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("input")
p.add_argument("root_cell")
p.add_argument("output")
p.add_argument("--groupby", default="leiden")
args = p.parse_args()
adata = sc.read_h5ad(args.input)
if "X_pca" not in adata.obsm or args.groupby not in adata.obs or adata.obs[args.groupby].nunique() < 2:
    p.error("input needs X_pca and at least two groups")
matches = np.flatnonzero(adata.obs_names == args.root_cell)
if len(matches) != 1:
    p.error("root barcode must identify exactly one cell")
sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca")
sc.tl.paga(adata, groups=args.groupby)
sc.tl.diffmap(adata, n_comps=15)
adata.uns["iroot"] = int(matches[0])
sc.tl.dpt(adata, n_dcs=10, n_branchings=0)
adata.write_h5ad(args.output)
