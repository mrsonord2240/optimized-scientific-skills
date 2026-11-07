# Partial charges

```python
import json
from rdkit import Chem
from scripts.gasteiger_charges import gasteiger_charges
smiles = open('input.smi').readline().strip()
mol = Chem.MolFromSmiles(smiles)
if mol is None or mol.GetNumAtoms() == 0: raise ValueError('Invalid molecule')
mol_h, charges = gasteiger_charges(mol)
json.dump({'smiles': smiles, 'atoms': [{'index': a.GetIdx(), 'element': a.GetSymbol(), 'charge': charges[a.GetIdx()]} for a in mol_h.GetAtoms()]}, open('charges.json', 'w'), indent=2)
```

The helper is `scripts/gasteiger_charges.py`; it adds hydrogens, rejects unsupported elements and checks finite charge balance. Atom charges follow the explicit-hydrogen molecule's atom order. Match the method to the downstream force field; see `references/charges.md`.

Done when charges, atom order and method limitations are recorded.
