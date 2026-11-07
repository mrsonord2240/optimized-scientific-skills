"""Run source MEBOCOST inference on an existing normalized AnnData."""
import argparse
import scanpy as sc
from mebocost import mebocost

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input")
    p.add_argument("config")
    p.add_argument("output_csv")
    p.add_argument("--groupby", default="cell_type")
    p.add_argument("--species", default="human")
    p.add_argument("--min-cells", type=int, default=50)
    args = p.parse_args()
    adata = sc.read_h5ad(args.input)
    if args.groupby not in adata.obs or adata.obs[args.groupby].isna().any():
        p.error("group metadata is missing")
    counts = adata.obs[args.groupby].value_counts()
    keep = counts[counts >= args.min_cells].index
    print("Excluded groups:", counts[counts < args.min_cells].to_dict())
    if len(keep) < 2:
        p.error("at least two groups must meet the minimum")
    adata = adata[adata.obs[args.groupby].isin(keep)].copy()
    model = mebocost.create_obj(adata=adata, group_col=args.groupby, condition_col=None,
        met_est="mebocost", config_path=args.config, species=args.species,
        cutoff_exp="auto", cutoff_met="auto", cutoff_prop=0.15, sensor_type="All", thread=8)
    result = model.infer_commu(n_shuffle=1000, seed=12345, Return=True,
        min_cell_number=10, pval_method="permutation_test_fdr", pval_cutoff=0.05, thread=None)
    result.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    main()
