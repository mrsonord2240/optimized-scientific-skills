from rdkit import Chem
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
