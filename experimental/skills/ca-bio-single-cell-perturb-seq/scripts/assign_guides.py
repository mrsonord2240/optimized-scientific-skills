"""Assign guides from an existing guide-count AnnData using source mixture method."""
import argparse
import scanpy as sc
import pertpy as pt
from scipy import sparse

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("input")
p.add_argument("output")
args = p.parse_args()
gdo = sc.read_h5ad(args.input)
if sparse.issparse(gdo.X):
    gdo.X = gdo.X.tocsr()
gdo.layers["counts"] = gdo.X.copy()
pt.pp.GuideAssignment().assign_mixture_model(gdo, assigned_guides_key="assigned_guide")
gdo.write_h5ad(args.output)
