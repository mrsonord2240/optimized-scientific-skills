'''Per-label marker agreement: do the cells given label L express L's marker genes as reference L cells do?

Why: the kNN/distance gate cannot see a query cell type that is MISSING from the reference but sits next
to a reference type (it is confidently mislabelled). Marker genes are independent of the embedding, so
this is the downstream check to run on every predicted label before trusting it.

Inputs: reference_labeled.h5ad (layers['counts'], obs['cell_type']), query_annotated.h5ad (layers['counts'],
obs['predicted_label'], written by scarches_annotation.py) and a REQUIRED --markers JSON {label: [genes]} from
prior knowledge. There is no reference-derived mode: markers derived from a reference that lacks the missing
neighbour type are generic and printed a false all-clear on the Monocytes hold-out.
Per label L: retained = mean over L's markers of (detection rate in query cells called L) / (detection rate
in reference L cells), capped at 1. A label is FLAGGED when retained < --min-retained (default 0.6).
Writes marker_check.tsv. Reference: scanpy 1.12, anndata 0.13 | Verify API if version differs
'''
import argparse, json
import numpy as np, pandas as pd, scanpy as sc

ap = argparse.ArgumentParser()
ap.add_argument('--reference', default='reference_labeled.h5ad')
ap.add_argument('--query', default='query_annotated.h5ad')
ap.add_argument('--markers', required=True, help='JSON {label: [genes]} from prior knowledge')
ap.add_argument('--min-retained', type=float, default=0.6)
ap.add_argument('--min-cells', type=int, default=20)
ap.add_argument('--out', default='marker_check.tsv')
a = ap.parse_args()

def prep(path):
    ad = sc.read_h5ad(path)
    ad.X = ad.layers['counts'].copy()
    sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    return ad

ref, qry = prep(a.reference), prep(a.query)
ref_labels = ref.obs['cell_type'].astype(str)
markers = json.load(open(a.markers))

def det(ad, mask, genes):
    x = ad[mask, genes].X
    x = x.toarray() if hasattr(x, 'toarray') else np.asarray(x)
    return (x > 0).mean(axis=0)

rows, skipped, n_unchecked = [], [], 0
pred = qry.obs['predicted_label'].astype(str)
for lab, n in pred.value_counts().items():
    if lab == 'Unknown':
        continue
    genes = [g for g in markers.get(lab, []) if g in qry.var_names and g in ref.var_names]
    if n < a.min_cells or not genes:
        skipped.append(f'{lab} (n={n}, {len(genes)} usable markers)'); n_unchecked += n
        continue
    d_ref = det(ref, (ref_labels == lab).values, genes)
    d_q = det(qry, (pred == lab).values, genes)
    keep = d_ref > 0
    if not keep.any():
        skipped.append(f'{lab} (markers never detected in reference)'); n_unchecked += n
        continue
    retained = float(np.minimum(1.0, d_q[keep] / d_ref[keep]).mean())
    rows.append({'label': lab, 'n_query': int(n), 'n_markers': int(keep.sum()), 'retained': round(retained, 3),
                 'flagged': retained < a.min_retained})
res = pd.DataFrame(rows)
res.to_csv(a.out, sep='\t', index=False)
print(res.to_string(index=False))
n_called = int((pred != 'Unknown').sum())
print('NOT CHECKED (no verdict, treat as unverified):', skipped or 'none')
print(f'Flagged labels: {res.loc[res.flagged, "label"].tolist() or "none"}  (retained < {a.min_retained}; '
      'a flagged label is untrustworthy even when the gate passed it)')
if n_unchecked:
    print(f'UNVERIFIED: {n_unchecked} of {n_called} labelled cells ({n_unchecked / n_called:.1%}) carry a label with no marker verdict. '
          'This is not a pass; report them as unverified.')
