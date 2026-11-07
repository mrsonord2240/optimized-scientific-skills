"""Compute source Mixscape signatures on existing annotated expression."""
import argparse
import scanpy as sc
import pertpy as pt
from scipy import sparse

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("input")
p.add_argument("output")
args = p.parse_args()
adata = sc.read_h5ad(args.input)
if not {"perturbation", "gene_target"}.issubset(adata.obs.columns):
    p.error("perturbation and gene_target metadata are required")
if sparse.issparse(adata.X):
    adata.X = adata.X.toarray()
ms = pt.tl.Mixscape()
ms.perturbation_signature(adata, pert_key="perturbation", control="NT", n_neighbors=20, batch_size=1000)
ms.mixscape(adata, pert_key="gene_target", control="NT", layer="X_pert")
adata.write_h5ad(args.output)
