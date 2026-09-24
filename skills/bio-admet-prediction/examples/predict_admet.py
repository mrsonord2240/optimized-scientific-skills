#!/usr/bin/env python3
'''
ADMET prediction and drug-likeness filtering.
'''
# Reference: rdkit 2024.03+, pandas 2.2+ | Verify API if version differs

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors, Lipinski
from rdkit.Chem.QED import qed
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


def load_admetlab_results(
    csv_path,
    *,
    smiles_column='SMILES',
    required_columns=None,
    require_uncertainty=True,
    task_id_column=None,
):
    '''Validate CSV output produced by the current official ADMETlab 3.0 workflow.

    The API owner publishes the current output contract. This loader makes the
    minimum safety checks locally: a structure column, any caller-required
    columns, at least one uncertainty-like column, and an optional task ID.
    '''
    csv_path = Path(csv_path)
    if not csv_path.is_file():
        raise FileNotFoundError(
            f'ADMETlab result file not found: {csv_path}. Submit through the '
            'current official API workflow before loading its output.'
        )
    try:
        results = pd.read_csv(csv_path)
    except pd.errors.EmptyDataError as exc:
        raise ValueError('ADMETlab result file contains no predictions') from exc
    if results.empty:
        raise ValueError('ADMETlab result file contains no predictions')

    expected = {smiles_column}
    if required_columns:
        expected.update(required_columns)
    if task_id_column:
        expected.add(task_id_column)
    missing = sorted(expected.difference(results.columns))
    if missing:
        raise ValueError(
            'ADMETlab result file is missing required columns: '
            + ', '.join(missing)
        )

    uncertainty_columns = [
        column for column in results.columns
        if any(token in column.lower() for token in
               ('uncertainty', 'confidence', 'interval', 'std', 'variance'))
    ]
    if require_uncertainty and not uncertainty_columns:
        raise ValueError(
            'ADMETlab result file has no uncertainty-like column. Check the '
            'current official output contract or pass require_uncertainty=False '
            'only for a documented computed-property-only export.'
        )
    return results


def predict_admet_ai(smiles, *, endpoints=None, model=None):
    '''Run ADMET-AI 2.x offline for valid SMILES after applicability gating.

    ``endpoints`` is optional because task names are package-version specific.
    It is deliberately validated rather than silently omitting misspelled tasks.
    Pass a constructed model to reuse weights across batches.
    '''
    if isinstance(smiles, str):
        smiles = [smiles]
    smiles = list(smiles)
    if not smiles:
        raise ValueError('Provide at least one SMILES string for ADMET-AI prediction')

    try:
        from admet_ai import ADMETModel
    except ImportError as exc:
        raise ImportError(
            'ADMET-AI 2.x is required for this route. Install admet-ai>=2,<3 '
            'and verify its installed task names.'
        ) from exc

    predictor = model or ADMETModel(include_physchem=False, drugbank_path=None)
    predictions = predictor.predict(smiles)
    if isinstance(predictions, dict):
        predictions = pd.DataFrame([predictions], index=smiles)

    if endpoints is None:
        return predictions
    missing = sorted(set(endpoints).difference(predictions.columns))
    if missing:
        raise ValueError(
            'Requested ADMET-AI endpoints are not available in this installed '
            'version: ' + ', '.join(missing) + '. Available examples: ' +
            ', '.join(map(str, predictions.columns[:12]))
        )
    return predictions.loc[:, list(endpoints)]


_ORGANIC_ATOMIC_NUMBERS = frozenset({1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53})


def _morgan_fingerprint(mol):
    return AllChem.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(mol)


def fit_similarity_threshold(reference_smiles, *, percentile=5):
    '''Derive a low-tail leave-one-out max-similarity threshold from a reference set.'''
    reference_mols = [Chem.MolFromSmiles(value) for value in reference_smiles]
    if len(reference_mols) < 2 or any(mol is None for mol in reference_mols):
        raise ValueError('Reference set needs at least two valid SMILES strings')
    if not 0 <= percentile <= 100:
        raise ValueError('percentile must be between 0 and 100')
    fingerprints = [_morgan_fingerprint(mol) for mol in reference_mols]
    leave_one_out = [
        max(DataStructs.TanimotoSimilarity(fp, other)
            for other_index, other in enumerate(fingerprints) if other_index != index)
        for index, fp in enumerate(fingerprints)
    ]
    return float(np.percentile(leave_one_out, percentile))


def assess_applicability_domain(
    smiles,
    reference_smiles,
    *,
    threshold,
    max_heavy_atoms=70,
):
    '''Gate small-molecule inputs and report reference-set maximum similarity.

    This is a transparent local similarity screen, not a claim about the model's
    true training domain. ``threshold`` must be derived for the relevant
    reference chemistry, for example with ``fit_similarity_threshold``.
    '''
    if not 0 <= threshold <= 1:
        raise ValueError('threshold must be a similarity value between 0 and 1')
    if max_heavy_atoms < 1:
        raise ValueError('max_heavy_atoms must be positive')
    if isinstance(smiles, str):
        smiles = [smiles]
    reference_mols = [Chem.MolFromSmiles(value) for value in reference_smiles]
    if not reference_mols or any(mol is None for mol in reference_mols):
        raise ValueError('Reference set must contain only valid SMILES strings')
    reference_fps = [_morgan_fingerprint(mol) for mol in reference_mols]

    rows = []
    for value in smiles:
        mol = Chem.MolFromSmiles(value)
        row = {'smiles': value, 'max_tanimoto': np.nan, 'decision': 'reject', 'reason': ''}
        if mol is None:
            row['reason'] = 'SMILES could not be parsed'
        elif len(Chem.GetMolFrags(mol)) != 1:
            row['reason'] = 'disconnected salts or mixtures are outside this route'
        elif any(atom.GetAtomicNum() not in _ORGANIC_ATOMIC_NUMBERS for atom in mol.GetAtoms()):
            row['reason'] = 'contains a metal or unsupported inorganic element'
        else:
            heavy_atoms = mol.GetNumHeavyAtoms()
            similarity = max(
                DataStructs.TanimotoSimilarity(_morgan_fingerprint(mol), reference_fp)
                for reference_fp in reference_fps
            )
            row['max_tanimoto'] = similarity
            if heavy_atoms > max_heavy_atoms:
                row['decision'] = 'manual_review'
                row['reason'] = f'{heavy_atoms} heavy atoms exceeds configured limit {max_heavy_atoms}'
            elif similarity < threshold:
                row['decision'] = 'manual_review'
                row['reason'] = f'max similarity {similarity:.3f} is below threshold {threshold:.3f}'
            else:
                row['decision'] = 'predict'
                row['reason'] = 'organic input passes local similarity screen'
        rows.append(row)
    return pd.DataFrame(rows)


def calculate_druglikeness(mol):
    '''Calculate drug-likeness properties.'''
    if mol is None:
        return None

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotatable = Lipinski.NumRotatableBonds(mol)

    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])

    return {
        'MW': mw,
        'LogP': logp,
        'HBD': hbd,
        'HBA': hba,
        'TPSA': tpsa,
        'RotatableBonds': rotatable,
        'QED': qed(mol),
        'LipinskiViolations': violations,
        'VeberCompliant': rotatable <= 10 and tpsa <= 140
    }


def flag_pains(molecules):
    '''
    Flag PAINS patterns without treating them as categorical exclusions.
    Returns (unflagged_molecules, flagged_with_descriptions).
    '''
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)

    clean = []
    flagged = []

    for mol in molecules:
        if mol is None:
            continue
        entry = catalog.GetFirstMatch(mol)
        if entry is None:
            clean.append(mol)
        else:
            flagged.append((mol, entry.GetDescription()))

    return clean, flagged


def flag_structural_alerts(molecules, catalogs=None):
    '''Flag structural alerts using multiple catalogs without discarding molecules.'''
    if catalogs is None:
        catalogs = [
            FilterCatalogParams.FilterCatalogs.PAINS,
            FilterCatalogParams.FilterCatalogs.BRENK,
            FilterCatalogParams.FilterCatalogs.NIH
        ]

    params = FilterCatalogParams()
    for cat in catalogs:
        params.AddCatalog(cat)
    catalog = FilterCatalog(params)

    clean = []
    alerts = []

    for mol in molecules:
        if mol is None:
            continue
        matches = catalog.GetMatches(mol)
        if not matches:
            clean.append(mol)
        else:
            alerts.append((mol, [m.GetDescription() for m in matches]))

    return clean, alerts


def annotate_compounds(molecules):
    '''Annotate compounds; project-specific ranking decisions belong downstream.'''
    results = []

    for mol in molecules:
        if mol is None:
            continue

        props = calculate_druglikeness(mol)
        if props is None:
            continue

        results.append((mol, props))

    return results


def batch_druglikeness(molecules):
    '''Calculate drug-likeness for multiple molecules.'''
    results = []
    for mol in molecules:
        props = calculate_druglikeness(mol)
        if props:
            props['SMILES'] = Chem.MolToSmiles(mol)
            results.append(props)
    return pd.DataFrame(results)


if __name__ == '__main__':
    smiles_list = ['CCO', 'c1ccccc1O', 'CC(=O)Oc1ccccc1C(=O)O']
    molecules = [Chem.MolFromSmiles(s) for s in smiles_list]

    print('Drug-likeness analysis:')
    df = batch_druglikeness(molecules)
    print(df[['SMILES', 'MW', 'LogP', 'QED', 'LipinskiViolations']].to_string())

    print('\nPAINS filtering:')
    unflagged, flagged = flag_pains(molecules)
    print(f'Unflagged: {len(unflagged)}, Flagged: {len(flagged)}')

    print('\nAnnotated compounds:')
    annotated = annotate_compounds(molecules)
    print(f'{len(annotated)} compounds annotated for project-specific review')
