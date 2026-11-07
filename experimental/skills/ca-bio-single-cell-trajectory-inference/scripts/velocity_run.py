"""Compute source deterministic velocity on existing splicing layers."""
import argparse
import scanpy as sc
import scvelo as scv

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input")
    p.add_argument("output")
    args = p.parse_args()
    adata = sc.read_h5ad(args.input)
    if not {"spliced", "unspliced"}.issubset(adata.layers):
        p.error("spliced and unspliced layers are required")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20)
    adata.layers["normalized_X"] = adata.X.copy()
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var["highly_variable"]].copy()
    adata.X = adata.layers.pop("normalized_X")
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(adata, mode="deterministic")
    scv.tl.velocity_graph(adata, n_jobs=1, show_progress_bar=False)
    scv.tl.velocity_confidence(adata)
    scv.tl.velocity_pseudotime(adata)
    adata.write_h5ad(args.output)

if __name__ == "__main__":
    main()
