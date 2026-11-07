"""Demultiplex public-format h5ad and aligned raw HTO count CSV."""
import argparse
import pandas as pd
import scanpy as sc
import scanpy.external as sce

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("gex")
p.add_argument("hto_csv")
p.add_argument("tags")
p.add_argument("output")
args = p.parse_args()
adata = sc.read_h5ad(args.gex)
hto = pd.read_csv(args.hto_csv, index_col=0)
tags = args.tags.split(",")
if len(tags) < 2 or len(set(tags)) != len(tags) or not set(tags).issubset(hto.columns):
    p.error("provide at least two unique HTO column names")
if not hto.index.is_unique or not adata.obs_names.isin(hto.index).all():
    p.error("HTO barcodes must be unique and cover GEX cells")
adata.obs[tags] = hto.loc[adata.obs_names, tags]
sce.pp.hashsolo(adata, cell_hashing_columns=tags, priors=(0.01, 0.8, 0.19),
                 number_of_noise_barcodes=1 if len(tags) <= 3 else None)
negative_fraction = (adata.obs["Classification"] == "Negative").mean()
if negative_fraction > 0.9:
    raise RuntimeError("Over 90% Negative; inspect noise, staining and priors")
adata.write_h5ad(args.output)
