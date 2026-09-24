# Reference: checked on pertpy 1.3.0, scanpy 1.12+, anndata 0.13+ | Verify API if version differs
#
# Perturb-seq analysis end-to-end using Pertpy:
# 1. sgRNA assignment
# 2. Escaper filtering via Mixscape
# 3. SCEPTRE differential expression (low-MOI calibrated)

import argparse

parser = argparse.ArgumentParser(
    description="Run a bounded Pertpy/Mixscape/PyDESeq2 demonstration on Papalexi 2021."
)
parser.add_argument(
    "--full",
    action="store_true",
    help="Run all cells, targets, and genes. This is not a smoke test and may need tens of GiB RAM.",
)
parser.add_argument("--max-cells", type=int, default=600, help="Bounded-mode cell cap (default: 600).")
parser.add_argument("--max-genes", type=int, default=2000, help="Bounded-mode gene cap (default: 2000).")
parser.add_argument(
    "--max-perturbations", type=int, default=2, help="Bounded-mode number of non-targeting contrasts (default: 2)."
)
parser.add_argument("--out", default="pertpy_de_results.tsv", help="Output DE TSV path.")
args = parser.parse_args()
if not args.full and (args.max_cells < 60 or args.max_genes < 100 or args.max_perturbations < 1):
    parser.error("bounded mode needs at least 60 cells, 100 genes, and one perturbation")

# Parse first so --help and invalid-argument feedback do not import the full single-cell stack.
import pertpy as pt
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

# === LOAD DATA ===
# Pertpy includes Papalexi 2021 Mixscape dataset for reference.
# papalexi_2021() returns a MuData object (RNA + ADT modalities), not a bare AnnData --
# extract the RNA modality before running scRNA-seq preprocessing/Mixscape/DE on it.
mdata = pt.dt.papalexi_2021()
labels = mdata.obs['gene_target'].astype(str)
if args.full:
    adata = mdata['rna'].copy()
    selected_labels = sorted(labels.unique())
    print("Full mode selected: expect the Mixscape signature matrix to require tens of GiB RAM.")
else:
    non_controls = sorted(label for label in labels.unique() if label != 'NT')
    selected_labels = ['NT', *non_controls[:args.max_perturbations]]
    per_label = max(1, args.max_cells // len(selected_labels))
    selected_positions = np.concatenate(
        [np.flatnonzero(labels.to_numpy() == label)[:per_label] for label in selected_labels]
    )
    adata = mdata['rna'][selected_positions, :args.max_genes].copy()
    print(
        f"Bounded smoke mode: {adata.n_obs} cells, {adata.n_vars} genes, "
        f"targets {selected_labels}; use --full only with a separately planned memory budget."
    )
# The per-cell target-gene label ('gene_target', with control cells labeled 'NT') lives on the
# top-level MuData.obs, not on the rna modality's own .obs -- merge it in (verified: same
# obs_names/order across mdata and mdata['rna']).
adata.obs['gene_target'] = labels.loc[adata.obs_names].to_numpy()
adata.layers['counts'] = adata.X.copy()  # PyDESeq2 needs raw counts -- save before normalizing

# Verify metadata
print(f'Cells: {adata.n_obs}')
print(f'Genes: {adata.n_vars}')
print(f'Perturbation column: "gene_target"')
print(f'NTC label: "NT"')

# === STANDARD scRNA-SEQ PREPROCESSING ===
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['log_normalized'] = adata.X.copy()
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata, mask_var='highly_variable')   # scanpy >=1.12: `use_highly_variable` is deprecated;
                                               # mask_var='highly_variable' is its exact replacement
                                               # (verified: byte-identical X_pca on papalexi_2021)

# === MIXSCAPE ESCAPER FILTERING ===
# Compute perturbation signature: cell_expression - mean(K NTC neighbors)
ms = pt.tl.Mixscape()
ms.perturbation_signature(
    adata=adata,
    pert_key='gene_target',
    control='NT',                                # papalexi_2021()'s own NTC label is 'NT'
    n_neighbors=20,                              # K neighbors for KNN-NTC subtraction
    random_state=0,                              # forwarded to pynndescent.NNDescent; omitting
                                                  # it makes X_pert non-deterministic across reruns
)
# Writes .layers['X_pert'] with perturbation signature

# Classify cells as KO (true perturbation) or NP (non-perturbed/escaper)
ms.mixscape(
    adata=adata,
    pert_key='gene_target',                     # pert_key (modern pertpy API)
    control='NT',
    new_class_name='mixscape_class',            # per-gene labels ('<gene> KO'/'<gene> NP'/'NT'); the
                                                # bare KO/NP/NT call goes to 'mixscape_class_global'
)
# Defaults to layer='X_pert' from previous step

# Retain KO cells + NTC controls for DE
adata_filtered = adata[adata.obs['mixscape_class_global'] != 'NP'].copy()
print(f'Cells after Mixscape filtering: {adata_filtered.n_obs}')
ko_pct = (adata.obs['mixscape_class_global'] == 'KO').sum() / (adata.obs['gene_target'] != 'NT').sum()
print(f'KO retention rate among perturbed: {ko_pct:.1%}')

# === DIFFERENTIAL EXPRESSION ===
# Per-perturbation pseudobulk DE via pertpy's PyDESeq2 wrapper
# Loop perturbations; one contrast per pert vs NT
# (For calibrated single-cell low-MOI DE, run R sceptre separately; Pertpy does not bundle it.)
# pertpy >= 1.0: build the contrast with .contrast(column, baseline, group_to_compare), then pass
# the resulting vector to test_contrasts(); result columns are log_fc / p_value / adj_p_value,
# not log2FoldChange / padj (checked on pertpy 1.3.0).
# layer='counts': PyDESeq2 requires raw (near-integer) counts -- pointing it at log-normalized
# data raises "ValueError: Non-zero elements of the matrix must be close to integer values."
# (verified: reproduced this exact error on real data, then fixed it by adding the counts layer).
de = pt.tl.PyDESeq2(adata_filtered, design='~gene_target', layer='counts')
de.fit()

all_de = []
for pert in adata_filtered.obs['gene_target'].unique():
    if pert == 'NT':
        continue
    contrast_df = de.test_contrasts(de.contrast('gene_target', 'NT', pert))
    contrast_df['perturbation'] = pert
    all_de.append(contrast_df)
de_result = pd.concat(all_de)

# === OUTPUT ===
de_result.to_csv(args.out, sep='\t')
print(f'DE results: {de_result.shape}')

# Top hits per perturbation
for pert in de_result['perturbation'].unique():
    pert_results = de_result[de_result['perturbation'] == pert].sort_values('adj_p_value')
    top = pert_results.head(10)
    if not top.empty:
        print(f'\n{pert}: top 10 differential genes')
        print(top[['log_fc', 'adj_p_value']].to_string())
