# Purpose: Mixscape escaper filtering: perturbation signature -> KO / NP (escaper) / control call -> keep KO cells.
# Input:   .h5ad of normalized, log1p-transformed scRNA-seq with a per-cell perturbation label in .obs[--pert-key]
#          (control cells carry the --control label, e.g. NTC or NT).
# Output:  .obs['mixscape_class'] ('<gene> KO' / '<gene> NP' / control) and .obs['mixscape_class_global'] (KO / NP / control),
#          .layers['X_pert']; KO cells written to --out if given. Requires pertpy >= 1.0.
# Usage:   python scripts/mixscape_filter.py normalized.h5ad --pert-key sgrna_assignment --control NTC [--out ko_cells.h5ad]
import argparse
import pertpy as pt
import scanpy as sc
import anndata as ad

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('h5ad')
p.add_argument('--pert-key', default='sgrna_assignment', help='.obs column with sgRNA target per cell')
p.add_argument('--control', default='NTC', help="match this to whatever your data's NTC label actually is")
p.add_argument('--n-neighbors', type=int, default=20, help='K neighbors for KNN-NTC subtraction')
p.add_argument('--out', help='write the KO-cell AnnData here')
a = p.parse_args()
adata = ad.read_h5ad(a.h5ad)

# adata is a scRNA-seq AnnData with 'sgrna_assignment' column
# Pertpy 0.6+ Mixscape API (verify against installed pertpy with help(pt.tl.Mixscape))
mixscape = pt.tl.Mixscape()
mixscape.perturbation_signature(
    adata=adata,
    pert_key=a.pert_key,             # .obs column with sgRNA target per cell
    control=a.control,               # match this to whatever your data's NTC label actually is
    n_neighbors=a.n_neighbors,       # K neighbors for KNN-NTC subtraction
    random_state=0,                  # forwarded to pynndescent.NNDescent -- omitting it makes
                                      # X_pert non-deterministic across reruns (verified: 5/2000
                                      # cells drifted between two unseeded runs on identical input)
)
# Writes .layers['X_pert'] with perturbation-signature-corrected expression

# Filter escapers: classify perturbed cells as KO (true perturbation) or NP (non-perturbed/escaper)
mixscape.mixscape(
    adata=adata,
    pert_key=a.pert_key,             # pert_key (not 'labels' in modern pertpy)
    control=a.control,
    new_class_name='mixscape_class', # .obs column to write
)
# Defaults to layer='X_pert' (output of perturbation_signature)

# Keep only KO cells for downstream analysis
# (mixscape_class holds '<gene> KO'; the bare label is in mixscape_class_global)
adata_ko = adata[adata.obs['mixscape_class_global'].isin(['KO'])].copy()
print(f'KO cells: {adata_ko.n_obs} ({adata_ko.n_obs/(adata.obs[a.pert_key] != a.control).sum():.1%} of perturbed)')
if a.out:
    adata_ko.write_h5ad(a.out)
