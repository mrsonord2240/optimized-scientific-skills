# Reference: RDKit 2024.09+, chembl_structure_pipeline 1.2+, pandas 2.2+ | Verify API if version differs

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from chembl_structure_pipeline import standardize_mol, get_parent_mol
import pandas as pd


def chembl_standardize(smi, multi_fragment_policy='flag'):
    '''Apply ChEMBL standardization and make ambiguous parents explicit.

    ``get_parent_mol`` can return a multi-fragment parent when all fragments are
    on ChEMBL's salt list.  The default flags that case instead of admitting an
    unstripped salt.  Set ``multi_fragment_policy='largest_fragment'`` only when
    the documented RDKit largest-organic-fragment fallback is appropriate.
    '''
    if multi_fragment_policy not in {'flag', 'largest_fragment'}:
        raise ValueError("multi_fragment_policy must be 'flag' or 'largest_fragment'")
    try:
        mol = Chem.MolFromSmiles(smi)
    except (TypeError, ValueError):
        return None, 'parse_failure'
    if mol is None:
        return None, 'parse_failure'
    try:
        std_mol = standardize_mol(mol)
        parent_mol, exclude = get_parent_mol(std_mol)
        if exclude:
            return None, 'excluded_by_chembl'
        if len(Chem.GetMolFrags(parent_mol)) > 1:
            if multi_fragment_policy == 'flag':
                return None, 'multi_fragment_parent'
            parent_mol = rdMolStandardize.LargestFragmentChooser(
                preferOrganic=True
            ).choose(parent_mol)
            return Chem.MolToSmiles(parent_mol), 'ok_largest_fragment_fallback'
        return Chem.MolToSmiles(parent_mol), 'ok'
    except Exception:
        return None, 'standardize_error'


def rdkit_standardize(smi, keep_isotopes=False):
    '''Full rdMolStandardize pipeline: sanitize -> largest fragment -> normalize -> uncharge -> canon tautomer.'''
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    Chem.SanitizeMol(mol)

    largest = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
    mol = largest.choose(mol)

    normalizer = rdMolStandardize.Normalizer()
    mol = normalizer.normalize(mol)

    # canonicalOrder=True selects neutralization sites in canonical order.
    uncharger = rdMolStandardize.Uncharger(canonicalOrder=True)
    mol = uncharger.uncharge(mol)

    enumerator = rdMolStandardize.TautomerEnumerator()
    mol = enumerator.Canonicalize(mol)

    if not keep_isotopes:
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return Chem.MolToSmiles(mol)


def prepare_qsar_data(
    df,
    smiles_col='smiles',
    activity_col='pIC50',
    pipeline='chembl',
    keep_isotopes_col=None,
    chembl_multi_fragment_policy='flag',
    max_activity_range=1.0,
    replicate_policy='flag',
):
    '''Prepare audited QSAR rows and return ``(deduplicated_rows, status_counts)``.

    ``keep_isotopes_col`` is an optional boolean column used per record.  A
    replicate group whose max-min activity span exceeds ``max_activity_range``
    is flagged by default; use ``replicate_policy='drop'`` only under a
    documented assay-QC policy.
    '''
    required = [smiles_col, activity_col]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required column(s): {', '.join(missing)}")
    if keep_isotopes_col is not None:
        if keep_isotopes_col not in df.columns:
            raise ValueError(f"Missing keep-isotopes column: {keep_isotopes_col}")
        if not pd.api.types.is_bool_dtype(df[keep_isotopes_col]):
            raise ValueError(f"{keep_isotopes_col} must have boolean dtype")
    if pipeline not in {'chembl', 'rdkit'}:
        raise ValueError("pipeline must be 'chembl' or 'rdkit'")
    if replicate_policy not in {'flag', 'drop'}:
        raise ValueError("replicate_policy must be 'flag' or 'drop'")
    if max_activity_range < 0:
        raise ValueError('max_activity_range must be non-negative')

    rows = []
    status_counts = {
        'parse_failure': 0,
        'excluded_by_chembl': 0,
        'multi_fragment_parent': 0,
        'ok_largest_fragment_fallback': 0,
        'standardize_error': 0,
        'inorganic_no_carbon': 0,
        'ok': 0,
    }
    for _, row in df.iterrows():
        smi = row[smiles_col]
        if pipeline == 'chembl':
            std_smi, status = chembl_standardize(
                smi, multi_fragment_policy=chembl_multi_fragment_policy
            )
        else:
            keep_isotopes = (
                bool(row[keep_isotopes_col]) if keep_isotopes_col is not None else False
            )
            try:
                std_smi = rdkit_standardize(smi, keep_isotopes=keep_isotopes)
            except Exception:
                std_smi = None
            status = 'ok' if std_smi else 'failed'
        if status == 'failed':
            status = 'standardize_error'
        if std_smi is None or not status.startswith('ok'):
            status_counts[status] = status_counts.get(status, 0) + 1
            continue
        mol = Chem.MolFromSmiles(std_smi)
        # Skip inorganic / no-carbon compounds
        if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) == 0:
            status_counts['inorganic_no_carbon'] += 1
            continue
        status_counts[status] = status_counts.get(status, 0) + 1
        rows.append({
            'smiles': std_smi,
            'inchikey': Chem.MolToInchiKey(mol),
            'activity': row[activity_col],
            'original_smiles': smi,
            'kept_isotopes': (
                bool(row[keep_isotopes_col]) if keep_isotopes_col is not None else False
            ),
        })
    df_std = pd.DataFrame(rows)
    if df_std.empty:
        empty = pd.DataFrame(columns=[
            'inchikey', 'smiles', 'activity', 'activity_std',
            'activity_range', 'replicate_disagreement', 'n_replicates',
            'original_smiles', 'kept_isotopes',
        ])
        return empty, status_counts

    # Deduplicate by InChIKey; preserve the observed disagreement for QC.
    df_dedup = df_std.groupby('inchikey').agg(
        smiles=('smiles', 'first'),
        activity=('activity', 'mean'),
        activity_std=('activity', 'std'),
        activity_range=('activity', lambda values: values.max() - values.min()),
        n_replicates=('activity', 'count'),
        original_smiles=('original_smiles', list),
        kept_isotopes=('kept_isotopes', 'any'),
    ).reset_index()
    df_dedup['replicate_disagreement'] = (
        (df_dedup['n_replicates'] > 1)
        & (df_dedup['activity_range'] > max_activity_range)
    )
    if replicate_policy == 'drop':
        df_dedup = df_dedup.loc[~df_dedup['replicate_disagreement']].reset_index(drop=True)
    return df_dedup, status_counts


if __name__ == '__main__':
    test_smiles = [
        'CC(=O)Oc1ccccc1C(=O)O',           # aspirin
        '[Na+].CC(=O)Oc1ccccc1C(=O)[O-]',  # aspirin sodium salt
        'CCC(=O)Nc1ccc(C(=O)C)cc1',         # different tautomer-ambiguous
        'invalid_smiles',                    # parse failure
    ]
    df = pd.DataFrame({'smiles': test_smiles, 'pIC50': [6.0, 6.2, 5.5, None]})
    result, status_counts = prepare_qsar_data(df)
    print('[status] ' + '  '.join(f'{key}={value}' for key, value in status_counts.items()))
    print(result[[
        'inchikey', 'smiles', 'activity', 'activity_range',
        'replicate_disagreement', 'n_replicates',
    ]])
