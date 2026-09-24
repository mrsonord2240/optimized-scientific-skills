# Reference: RDKit 2024.09+, pandas 2.2+, scikit-learn 1.4+ | Verify API if version differs

from collections import defaultdict
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import pandas as pd
import random


def get_bemis_murcko_scaffold(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold)


def get_generic_framework(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    framework = MurckoScaffold.MakeScaffoldGeneric(scaffold)
    return Chem.MolToSmiles(framework)


def scaffold_clusters(df, smiles_col='smiles'):
    '''Group zero-based row positions by Bemis-Murcko scaffold.'''
    clusters = defaultdict(list)
    invalid_positions = []
    for pos, smi in enumerate(df[smiles_col].tolist()):
        scaff = get_bemis_murcko_scaffold(smi)
        if scaff is None:
            invalid_positions.append(pos)
        else:
            clusters[scaff].append(pos)
    if invalid_positions:
        raise ValueError(f'Invalid SMILES at row positions: {invalid_positions}')
    return dict(clusters)


def scaffold_split(df, smiles_col='smiles', train_frac=0.8, seed=42,
                   return_diagnostics=False):
    '''Assign whole scaffolds and expose split diagnostics when requested.'''
    clusters = scaffold_clusters(df, smiles_col)
    if not 0 < train_frac < 1:
        raise ValueError('train_frac must be strictly between 0 and 1')
    scaffold_sets = list(clusters.values())
    n_total = sum(len(s) for s in scaffold_sets)
    n_train = int(n_total * train_frac)
    n_test = n_total - n_train

    rng = random.Random(seed)
    rng.shuffle(scaffold_sets)
    scaffold_sets.sort(key=lambda group: len(group), reverse=True)

    if len(scaffold_sets) < 2:
        raise ValueError('A scaffold split requires at least two distinct scaffolds')
    if n_train == 0 or n_test == 0:
        raise ValueError('train_frac produces an empty partition; use a larger dataset or a less extreme fraction')

    train_idx, test_idx = [], []
    train_group_sizes, test_group_sizes = [], []
    for scaff_set in scaffold_sets:
        if len(train_idx) / n_train <= len(test_idx) / n_test:
            train_idx.extend(scaff_set)
            train_group_sizes.append(len(scaff_set))
        else:
            test_idx.extend(scaff_set)
            test_group_sizes.append(len(scaff_set))

    # Exact fractions are possible only when singleton groups can close the gap.
    train_singletons = [s for s in scaffold_sets if len(s) == 1 and s[0] in train_idx]
    test_singletons = [s for s in scaffold_sets if len(s) == 1 and s[0] in test_idx]
    while len(train_idx) < n_train and test_singletons:
        group = test_singletons.pop()
        test_idx.remove(group[0]); train_idx.extend(group)
        test_group_sizes.remove(1); train_group_sizes.append(1)
    while len(train_idx) > n_train and train_singletons:
        group = train_singletons.pop()
        train_idx.remove(group[0]); test_idx.extend(group)
        train_group_sizes.remove(1); test_group_sizes.append(1)

    diagnostics = {
        'requested_train_fraction': train_frac,
        'achieved_train_fraction': len(train_idx) / n_total,
        'scaffold_overlap': 0,
        'train_singleton_compound_fraction': sum(size == 1 for size in train_group_sizes) / len(train_idx),
        'test_singleton_compound_fraction': sum(size == 1 for size in test_group_sizes) / len(test_idx),
    }
    if any(len(group) > 1 for group in scaffold_sets):
        assert diagnostics['test_singleton_compound_fraction'] < 1, (
            'All test compounds are singleton scaffolds; inspect the split diagnostics')
    train = df.iloc[train_idx].reset_index(drop=True)
    test = df.iloc[test_idx].reset_index(drop=True)
    if return_diagnostics:
        return train, test, diagnostics
    return train, test


def detect_analog_series(df, smiles_col='smiles', min_size=3):
    '''Identify series; three members is a repository starting definition.'''
    clusters = scaffold_clusters(df, smiles_col)
    series = {scaff: cmpds for scaff, cmpds in clusters.items()
              if len(cmpds) >= min_size}
    summary = pd.DataFrame([
        {'scaffold': scaff, 'n_members': len(cmpds)}
        for scaff, cmpds in series.items()
    ]).sort_values('n_members', ascending=False)
    return series, summary


if __name__ == '__main__':
    sample_smiles = [
        'c1ccc(C(=O)NCC)cc1F',
        'c1ccc(C(=O)NCCC)cc1Cl',
        'c1ccc(C(=O)NCC)cc1Br',
        'CCC(=O)N1CCN(c2ccc(F)cc2)CC1',
        'CCC(=O)N1CCN(c2ccc(Cl)cc2)CC1',
    ]
    df = pd.DataFrame({'smiles': sample_smiles, 'pIC50': [6.5, 6.2, 6.0, 7.1, 7.0]})
    # Small-demo values expose both partitions reproducibly; production ratios
    # and minimum series size must follow deployment and dataset size.
    train, test = scaffold_split(df, train_frac=0.6, seed=42)
    print(f'Train: {len(train)}, Test: {len(test)}')
    series, summary = detect_analog_series(df, min_size=2)
    print(summary.to_string())
