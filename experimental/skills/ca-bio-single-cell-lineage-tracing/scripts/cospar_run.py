"""Infer a CoSpar transition map from an existing clonal AnnData object."""
import argparse
import scanpy as sc
import cospar as cs

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("input")
p.add_argument("output")
p.add_argument("dataset_label")
args = p.parse_args()
adata = sc.read_h5ad(args.input)
if "X_clone" not in adata.obsm or "time_info" not in adata.obs:
    p.error("input requires X_clone and time_info")
emb = adata.copy()
sc.pp.normalize_total(emb)
sc.pp.log1p(emb)
sc.pp.pca(emb)
sc.pp.neighbors(emb)
sc.tl.umap(emb)
adata = cs.pp.initialize_adata_object(adata, X_clone=adata.obsm["X_clone"],
    time_info=adata.obs["time_info"], X_pca=emb.obsm["X_pca"],
    X_emb=emb.obsm["X_umap"], data_des=args.dataset_label)
adata = cs.tmap.infer_Tmap_from_multitime_clones(
    adata, smooth_array=[15, 10, 5], sparsity_threshold=0.1)
adata.write_h5ad(args.output)
