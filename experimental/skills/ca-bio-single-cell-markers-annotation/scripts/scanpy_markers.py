"""Rank cluster markers descriptively using Scanpy."""
import argparse
import scanpy as sc

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("input")
p.add_argument("output_csv")
p.add_argument("--groupby", default="leiden")
args = p.parse_args()
adata = sc.read_h5ad(args.input)
if args.groupby not in adata.obs or adata.obs[args.groupby].nunique() < 2:
    p.error("groupby must contain at least two clusters")
sc.tl.rank_genes_groups(adata, groupby=args.groupby, method="wilcoxon", pts=True,
                        corr_method="benjamini-hochberg", use_raw=False)
sc.get.rank_genes_groups_df(adata, group=None).to_csv(args.output_csv, index=False)
