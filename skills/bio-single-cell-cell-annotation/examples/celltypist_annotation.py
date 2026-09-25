# Reference: pandas 2.2+, scanpy 1.10+, scikit-learn 1.4+ | Verify API if version differs
import argparse
from pathlib import Path

import scanpy as sc
import celltypist
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Annotate seeded single-cell clusters with a cached CellTypist model.')
parser.add_argument('--input', default='clustered.h5ad', help='AnnData with raw counts, seeded leiden labels, and UMAP')
parser.add_argument('--model', default='Immune_All_Low.pkl', help='Locally cached CellTypist model name or path')
parser.add_argument('--output', default='adata_annotated.h5ad')
parser.add_argument('--figure', default='celltypist_annotation.png')
parser.add_argument('--counts-output', default='cell_type_counts.csv')
args = parser.parse_args()

adata = sc.read_h5ad(args.input)
if 'counts' not in adata.layers:
    raise ValueError("Input must retain raw counts in adata.layers['counts'] for CP10K-log1p normalization.")
if 'leiden' not in adata.obs:
    raise ValueError("Input must contain an intentionally seeded over-clustering in adata.obs['leiden'].")
if 'X_umap' not in adata.obsm:
    raise ValueError("Input must contain adata.obsm['X_umap'] for the annotation plots.")

adata.X = adata.layers['counts'].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Resolve the model without CellTypist's model-index helpers: Model.load() may
# consult and populate the default cache when it receives a bare model name.
# Passing an existing explicit path keeps acquisition separate from execution.
model_path = Path(args.model).expanduser()
if not model_path.is_file():
    model_path = Path(celltypist.models.models_path) / args.model
if not model_path.is_file():
    raise FileNotFoundError(
        f"CellTypist model is not cached: {args.model}. Download it separately or pass an existing file with --model."
    )
model = celltypist.models.Model.load(model=model_path.resolve().as_posix())
predictions = celltypist.annotate(
    adata,
    model=model,
    majority_voting=True,
    over_clustering='leiden',
)
adata = predictions.to_adata()

adata.obs['cell_type'] = adata.obs['majority_voting']
adata.obs['annotation_confidence'] = adata.obs['conf_score']

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color='cell_type', ax=axes[0], show=False, title='Cell Types')
sc.pl.umap(adata, color='annotation_confidence', ax=axes[1], show=False,
           title='Confidence Score', cmap='viridis')
plt.tight_layout()
plt.savefig(args.figure, dpi=150, bbox_inches='tight')
plt.close()

confidence_threshold = 0.5
adata.obs['high_confidence'] = adata.obs['annotation_confidence'] > confidence_threshold
low_conf_cells = adata.obs[~adata.obs['high_confidence']]
print(f'Low confidence cells: {len(low_conf_cells)} ({len(low_conf_cells)/len(adata)*100:.1f}%)')

cell_type_counts = adata.obs['cell_type'].value_counts()
cell_type_counts.to_csv(args.counts_output)

adata.write(args.output)
