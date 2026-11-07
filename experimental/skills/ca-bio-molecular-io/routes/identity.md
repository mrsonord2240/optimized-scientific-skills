# InChI identity

```python
from rdkit import Chem
from rdkit.Chem.inchi import MolToInchi
mol = Chem.MolFromSmiles(open('input.smi').readline().strip())
if mol is None or mol.GetNumAtoms() == 0: raise ValueError('Invalid molecule')
open('identity.txt', 'w').write(MolToInchi(mol) + '\n')
```

Standard InChI can normalize mobile-hydrogen tautomers; for FixedH and collision investigation use `references/conversion.md`.

Done when identifiers and standardization settings are saved.
