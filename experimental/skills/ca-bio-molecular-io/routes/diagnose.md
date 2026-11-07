# Representation diagnosis

```python
from rdkit import Chem
mol = Chem.MolFromSmiles(open('input.smi').readline().strip(), sanitize=False)
if mol is None or mol.GetNumAtoms() == 0: raise ValueError('Invalid or empty SMILES grammar')
print(Chem.SanitizeMol(mol, catchErrors=True))
```

Separate syntax from sanitization failures. For aromaticity, stereo, V3000, PDB bond-order templates and charge-model caveats read `references/formats.md` and `references/failures.md`.

Done when failures and representation differences are identified.
