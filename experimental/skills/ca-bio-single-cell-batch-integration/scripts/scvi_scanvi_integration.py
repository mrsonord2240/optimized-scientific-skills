#!/usr/bin/env python
"""Seeded scVI/scANVI integration. Usage: scvi_scanvi_integration.py INPUT.h5ad OUTPUT.h5ad --labels-key cell_type."""
import argparse
import scanpy as sc
import scvi


def main():
    parser = argparse.ArgumentParser(description="Integrate raw counts with scVI and optional scANVI.")
    parser.add_argument("input_h5ad")
    parser.add_argument("output_h5ad")
    parser.add_argument("--batch-key", default="batch")
    parser.add_argument("--labels-key")
    parser.add_argument("--max-epochs", type=int, default=100)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    adata = sc.read_h5ad(args.input_h5ad)
    if args.batch_key not in adata.obs:
        raise KeyError(f"Missing batch key: {args.batch_key}")
    if args.labels_key and args.labels_key not in adata.obs:
        raise KeyError(f"Missing labels key: {args.labels_key}")
    scvi.settings.seed = args.seed
    adata.layers["counts"] = adata.X.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",
                                layer="counts", batch_key=args.batch_key)
    adata = adata[:, adata.var.highly_variable].copy()
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=args.batch_key)
    model = scvi.model.SCVI(adata, n_latent=10, gene_likelihood="zinb")
    model.train(max_epochs=args.max_epochs, accelerator="cpu")
    adata.obsm["X_scVI"] = model.get_latent_representation()
    if args.labels_key:
        labels = adata.obs[args.labels_key].astype(str)
        if "Unknown" not in set(labels):
            raise ValueError(f"{args.labels_key} must contain 'Unknown' to run scANVI")
        scanvi = scvi.model.SCANVI.from_scvi_model(model, "Unknown", labels_key=args.labels_key)
        scanvi.train(max_epochs=args.max_epochs, accelerator="cpu")
        adata.obsm["X_scANVI"] = scanvi.get_latent_representation()
        adata.obs["scanvi_label"] = scanvi.predict()
    adata.write_h5ad(args.output_h5ad)


if __name__ == "__main__":
    main()
