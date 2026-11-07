#!/usr/bin/env python
"""Seeded Harmony integration. Usage: harmony_integration.py INPUT.h5ad OUTPUT.h5ad [--batch-key batch] [--seed 0]."""
import argparse
import numpy as np
import scanpy as sc
import scanpy.external as sce


def main():
    parser = argparse.ArgumentParser(description="Integrate batches with Harmony in PCA space.")
    parser.add_argument("input_h5ad")
    parser.add_argument("output_h5ad")
    parser.add_argument("--batch-key", default="batch")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    adata = sc.read_h5ad(args.input_h5ad)
    if args.batch_key not in adata.obs:
        raise KeyError(f"Missing batch key: {args.batch_key}")
    np.random.seed(args.seed)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=args.batch_key)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, random_state=args.seed)
    sce.pp.harmony_integrate(adata, key=args.batch_key, random_state=args.seed)
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", random_state=args.seed)
    sc.tl.umap(adata, random_state=args.seed)
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False, random_state=args.seed)
    adata.write_h5ad(args.output_h5ad)


if __name__ == "__main__":
    main()
