"""Apply an explicitly reviewed cluster-to-label mapping."""
import argparse
import pandas as pd
import scanpy as sc

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("input")
p.add_argument("labels_csv")
p.add_argument("output")
p.add_argument("--groupby", default="leiden")
args = p.parse_args()
adata = sc.read_h5ad(args.input)
labels = pd.read_csv(args.labels_csv, dtype=str)
if not {"cluster", "label"}.issubset(labels.columns) or labels["cluster"].duplicated().any() or labels[["cluster", "label"]].isna().any().any():
    p.error("labels need unique, nonmissing cluster and label columns")
if args.groupby not in adata.obs:
    p.error("groupby column missing")
mapping = dict(zip(labels["cluster"], labels["label"]))
adata.obs["cell_type"] = adata.obs[args.groupby].astype(str).map(mapping).fillna("Unassigned")
adata.write_h5ad(args.output)
