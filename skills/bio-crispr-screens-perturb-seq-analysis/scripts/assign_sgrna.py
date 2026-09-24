# Purpose: per-cell sgRNA assignment by read-count threshold (single sgRNA / 'multiplet' / 'none').
# Input:   .h5ad whose .layers[<layer>] is a cells x sgRNAs count matrix and whose var_names are sgRNA names.
# Output:  the AnnData with .obs["sgrna_assignment"] added (written to --out if given); prints the assignment fractions.
# Usage:   python scripts/assign_sgrna.py sgrnas.h5ad [--layer sgrna_counts] [--threshold 10] [--out assigned.h5ad]
#          or import: from assign_sgrna import assign_sgrna
import argparse


# sgRNA assignment via threshold counting
def assign_sgrna(adata, sgrna_counts_layer='sgrna_counts', threshold=10):
    '''Per-cell sgRNA assignment. Returns single assignment or 'multiplet'/'none'.'''
    import numpy as np
    counts = adata.layers[sgrna_counts_layer]  # cells x sgRNAs
    # Works for a dense array or a scipy sparse layer (the usual 10X / FeatureBarcode case).
    # sparse sum/argmax return an (cells, 1) matrix, so ravel both -- a 2-D index array raises
    # pyarrow.lib.ArrowInvalid when it reaches .obs, and a 2-D mask raises AxisError from .sum().
    n_sgrna_per_cell = np.asarray((counts >= threshold).sum(axis=1)).ravel()
    top_sgrna = np.asarray(adata.var_names)[np.asarray(counts.argmax(axis=1)).ravel()]
    assignments = np.where(
        n_sgrna_per_cell == 0, 'none',
        np.where(n_sgrna_per_cell == 1, top_sgrna, 'multiplet'))
    adata.obs['sgrna_assignment'] = assignments
    return adata


if __name__ == '__main__':
    import anndata as ad
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('h5ad', help='AnnData with a cells x sgRNAs count layer')
    p.add_argument('--layer', default='sgrna_counts')
    p.add_argument('--threshold', type=int, default=10)
    p.add_argument('--out', help='write the annotated AnnData here')
    a = p.parse_args()
    adata = assign_sgrna(ad.read_h5ad(a.h5ad), a.layer, a.threshold)
    print(adata.obs['sgrna_assignment'].value_counts(normalize=True).round(3).to_string())
    if a.out:
        adata.write_h5ad(a.out)
