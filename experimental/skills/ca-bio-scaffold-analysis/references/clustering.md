## Library Chemotype Clustering

**Goal:** Group compounds by shared Bemis-Murcko scaffold.

**Approach:** Compute scaffold for each compound; group by scaffold SMILES.

```python
from collections import defaultdict

def scaffold_clusters(smiles_list):
    clusters = defaultdict(list)
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_smi = Chem.MolToSmiles(scaffold)
        clusters[scaffold_smi].append(smi)
    return clusters
```

Output: dict {scaffold_smiles: [compound_smiles, ...]}. Cluster sizes inform library diversity.

## Bemis-Murcko Scaffold Split (ML)

For QSAR / ML, random train/test split causes data leakage: compounds from the same chemotype (analogs in same series) end up in both. Bemis-Murcko split puts entire scaffolds in train or test, never both.

```python
from rdkit.Chem.Scaffolds import MurckoScaffold

def scaffold_split(df, smiles_col='smiles', train_frac=0.8, seed=42,
                   return_diagnostics=False):
    import random
    rng = random.Random(seed)

    scaffolds = defaultdict(list)
    invalid_positions = []
    for pos, smi in enumerate(df[smiles_col].tolist()):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            invalid_positions.append(pos)
            continue
        scaff = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
        scaffolds[scaff].append(pos)

    if invalid_positions:
        raise ValueError(f'Invalid SMILES at row positions: {invalid_positions}')

    if not 0 < train_frac < 1:
        raise ValueError('train_frac must be strictly between 0 and 1')
    scaffold_sets = list(scaffolds.values())
    # Shuffle before size sorting so equal-size groups do not get a positional bias.
    rng.shuffle(scaffold_sets)
    scaffold_sets.sort(key=lambda x: len(x), reverse=True)

    n_total = sum(len(s) for s in scaffold_sets)
    n_train = int(n_total * train_frac)
    if len(scaffold_sets) < 2:
        raise ValueError('A scaffold split requires at least two scaffolds')
    n_test = n_total - n_train
    if n_train == 0 or n_test == 0:
        raise ValueError('train_frac produces an empty partition; use a larger dataset or a less extreme fraction')
    train_idx, test_idx = [], []
    train_group_sizes, test_group_sizes = [], []
    # Allocate large groups to the currently less-filled partition. This preserves
    # whole-scaffold isolation while spreading multi-member groups across splits;
    # singleton groups then make the requested fraction as close as possible.
    for scaff_set in scaffold_sets:
        train_fill = len(train_idx) / n_train
        test_fill = len(test_idx) / n_test
        if train_fill <= test_fill:
            train_idx.extend(scaff_set)
            train_group_sizes.append(len(scaff_set))
        else:
            test_idx.extend(scaff_set)
            test_group_sizes.append(len(scaff_set))

    # A last singleton-only adjustment can hit an exact target without moving a
    # multi-member scaffold. Otherwise report the closest whole-group split.
    if len(train_idx) != n_train:
        # Move singleton rows only; every singleton is its own scaffold group.
        train_singleton_groups = [s for s in scaffold_sets if len(s) == 1 and s[0] in train_idx]
        test_singleton_groups = [s for s in scaffold_sets if len(s) == 1 and s[0] in test_idx]
        while len(train_idx) < n_train and test_singleton_groups:
            group = test_singleton_groups.pop()
            test_idx.remove(group[0]); train_idx.extend(group)
            test_group_sizes.remove(1); train_group_sizes.append(1)
        while len(train_idx) > n_train and train_singleton_groups:
            group = train_singleton_groups.pop()
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
    train, test = df.iloc[train_idx], df.iloc[test_idx]
    if return_diagnostics:
        return train, test, diagnostics
    return train, test
```

**Effect on benchmark metrics:** A scaffold split often produces different performance from a random split because it tests transfer across scaffold groups. The size and meaning of the gap are dataset- and deployment-dependent; it is not a direct universal measure of memorization.

**Caveat:** Bemis-Murcko split is *one* scaffold-split; for production ML, consider time split (newer compounds in test) or activity-cliff-balanced split.

**Class-imbalanced datasets:** The bundled splitter spreads large scaffold groups across the two partitions and returns `test_singleton_compound_fraction` with `return_diagnostics=True`; audit it with scaffold overlap and endpoint balance. Chemprop's `scaffold_balanced` split balances scaffold-group sizes; it is not label-stratified. If both group isolation and label balance are required, use a validated group-aware stratification procedure such as `StratifiedGroupKFold` where its assumptions fit, then audit every fold for scaffold overlap and endpoint balance.
