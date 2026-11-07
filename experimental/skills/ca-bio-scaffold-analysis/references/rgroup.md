## R-Group Decomposition

**Goal:** Given a defined scaffold and a set of analog compounds, extract the R-group at each numbered attachment point into a tabular SAR matrix.

```python
from rdkit.Chem import rdRGroupDecomposition as rgd

def decompose_series(compounds, scaffold_smiles_with_R):
    scaffold = Chem.MolFromSmiles(scaffold_smiles_with_R)
    if scaffold is None:
        raise ValueError('Invalid scaffold SMARTS/SMILES')
    parsed = [(i, Chem.MolFromSmiles(s)) for i, s in enumerate(compounds)]
    invalid = [i for i, mol in parsed if mol is None]
    if invalid:
        raise ValueError(f'Invalid compound SMILES at positions: {invalid}')
    mols = [mol for _, mol in parsed]
    decomp, unmatched = rgd.RGroupDecompose([scaffold], mols, asSmiles=True)
    unmatched_set = set(unmatched)
    matched_positions = [i for i in range(len(mols)) if i not in unmatched_set]
    return decomp, matched_positions, list(unmatched)

scaffold = 'c1ccc(C(=O)N[*:1])cc1-[*:2]'
compounds = ['c1ccc(C(=O)NCC)cc1F', 'c1ccc(C(=O)NCCC)cc1Cl']
table = decompose_series(compounds, scaffold)
```

Output: list of {'Core': scaffold, 'R1': r1_smiles, 'R2': r2_smiles} dicts. Used for Free-Wilson analysis (see reaction-enumeration skill).
