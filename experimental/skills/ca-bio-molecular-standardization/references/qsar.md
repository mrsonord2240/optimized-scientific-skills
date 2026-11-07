## Standardization for ML Training (avoiding data leakage)

**Goal:** Build a standardized + deduplicated training set with replicate-averaged activity for QSAR or ADMET model training.

**Approach:** Standardize every SMILES through the ChEMBL pipeline, compute InChIKey as canonical identity, group by InChIKey, and mean-aggregate activities; report replicate count for confidence weighting.

```python
import pandas as pd
from rdkit import Chem
from chembl_structure_pipeline import standardize_mol, get_parent_mol

def prepare_qsar_data(
    df, smiles_col='smiles', activity_col='pIC50',
    max_activity_range=1.0, replicate_policy='flag',
):
    missing = [name for name in (smiles_col, activity_col) if name not in df.columns]
    if missing:
        raise ValueError(f"Missing required column(s): {', '.join(missing)}")
    if replicate_policy not in {'flag', 'drop'}:
        raise ValueError("replicate_policy must be 'flag' or 'drop'")

    standardized = []
    status_counts = {
        'parse_failure': 0, 'excluded_by_chembl': 0,
        'multi_fragment_parent': 0, 'standardize_error': 0,
        'inorganic_no_carbon': 0, 'ok': 0,
    }
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(row[smiles_col])
        if mol is None:
            status_counts['parse_failure'] += 1
            continue
        try:
            mol = standardize_mol(mol)
            mol, exclude = get_parent_mol(mol)
            if exclude:
                status_counts['excluded_by_chembl'] += 1
                continue
            if len(Chem.GetMolFrags(mol)) > 1:
                status_counts['multi_fragment_parent'] += 1
                continue
        except Exception:
            status_counts['standardize_error'] += 1
            continue
        if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
            status_counts['inorganic_no_carbon'] += 1
            continue
        status_counts['ok'] += 1
        standardized.append({
            'smiles': Chem.MolToSmiles(mol),
            'inchikey': Chem.MolToInchiKey(mol),
            'activity': row[activity_col],
        })

    df_std = pd.DataFrame(standardized)
    if df_std.empty:
        return pd.DataFrame(), status_counts
    df_std = df_std.groupby('inchikey').agg(
        smiles=('smiles', 'first'),
        activity=('activity', 'mean'),
        activity_range=('activity', lambda values: values.max() - values.min()),
        n_replicates=('activity', 'count'),
    ).reset_index()
    df_std['replicate_disagreement'] = (
        (df_std['n_replicates'] > 1)
        & (df_std['activity_range'] > max_activity_range)
    )
    if replicate_policy == 'drop':
        df_std = df_std.loc[~df_std['replicate_disagreement']].reset_index(drop=True)
    return df_std, status_counts
```

The bundled `scripts/standardize_library.py` is the complete reference implementation. It additionally accepts `keep_isotopes_col`, a boolean per-record flag for mixed tracer libraries, and `chembl_multi_fragment_policy='largest_fragment'` when that fallback has been authorized.

`prepare_qsar_data` returns `(deduplicated_dataframe, status_counts)`. Emit and retain `status_counts`; it makes each excluded input traceable rather than silently dropping rows. Standard InChIKey may collapse some mobile-hydrogen tautomer representations, but this is not a substitute for an explicitly chosen tautomer policy.

`activity_range` is the max-minus-min value across every identity group's reported activities. The default threshold, `max_activity_range=1.0` log unit, sets `replicate_disagreement=True`; this is a review flag, not a reason to erase data automatically. Set `replicate_policy='drop'` only after documenting assay format, endpoint, units, and the agreed QC rule. ChEMBL exports can mix binding and functional assays, so reconcile assay scope before averaging.
