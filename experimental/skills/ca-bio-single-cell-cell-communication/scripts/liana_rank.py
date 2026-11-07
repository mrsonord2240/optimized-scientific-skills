"""Run source LIANA consensus ranking on an annotated h5ad."""
import argparse
import scanpy as sc
import liana as li

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("input")
p.add_argument("output_csv")
p.add_argument("--groupby", default="cell_type")
p.add_argument("--n-perms", type=int, default=1000)
args = p.parse_args()
adata = sc.read_h5ad(args.input)
if args.groupby not in adata.obs or adata.obs[args.groupby].nunique() < 2:
    p.error("groupby must name at least two cell types")
li.mt.rank_aggregate(adata, groupby=args.groupby, resource_name="consensus",
                     expr_prop=0.1, use_raw=False, n_perms=args.n_perms, verbose=True)
adata.uns["liana_res"].to_csv(args.output_csv, index=False)
