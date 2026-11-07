## Reading SMILES with Stereo Preservation

**Goal:** Parse SMILES while preserving stereo and aromatic-flag consistency.

**Approach:** Use `Chem.MolFromSmiles(smi)` with sanitization on, verify with round-trip canonicalization, and set explicit stereochemistry where the toolkit's perception missed it.

```python
from rdkit import Chem
from rdkit.Chem import AllChem

def parse_smiles_safe(smi):
    if not isinstance(smi, str) or not smi.strip():
        return None, 'empty_input'
    mol = Chem.MolFromSmiles(smi)
    if mol is None or mol.GetNumAtoms() == 0:
        return None, 'parse_failure'
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    canon = Chem.MolToSmiles(mol)
    round_trip = Chem.MolFromSmiles(canon)
    if Chem.MolToSmiles(round_trip) != canon:
        return mol, 'round_trip_unstable'
    return mol, 'ok'
```

## Reading SDF with Property Carryover

**Goal:** Load a multi-record SDF preserving per-molecule properties (Name, ID, IC50, etc.) used by downstream filtering and ML labeling.

**Approach:** Iterate via `SDMolSupplier(removeHs=False, sanitize=True)`, filter `None` (parse failures), and capture properties via `mol.GetPropsAsDict()`.

```python
from rdkit import Chem

supplier = Chem.SDMolSupplier('library.sdf', removeHs=False, sanitize=True)
mols = []
fails = []
for i, mol in enumerate(supplier):
    if mol is None or mol.GetNumAtoms() == 0:
        fails.append(i)
        continue
    props = mol.GetPropsAsDict()
    mols.append((mol, props))
print(f'parsed: {len(mols)}; failed: {len(fails)}')
```

An empty or whitespace-only SMILES is not a molecule: RDKit can return a zero-atom `Mol` rather than `None`, so reject it before parsing and check `GetNumAtoms() == 0` at file boundaries.

For a failed non-empty SMILES, first distinguish grammar from chemistry. Bad parentheses, unclosed rings, and invalid tokens return `None` even with `sanitize=False`; correct the source text. If `Chem.MolFromSmiles(smi, sanitize=False)` returns a molecule, use `Chem.SanitizeMol(mol, catchErrors=True)` to identify a chemistry failure such as kekulization, valence, or aromaticity.
