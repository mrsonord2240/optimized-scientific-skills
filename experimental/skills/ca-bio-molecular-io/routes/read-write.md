# Read and write structures

```python
from rdkit import Chem
from scripts.molecule_io import write_smiles
supplier = Chem.SDMolSupplier('input.sdf', removeHs=False, sanitize=True)
valid, failed = [], []
for i, mol in enumerate(supplier):
    (failed if mol is None or mol.GetNumAtoms() == 0 else valid).append(i if mol is None or mol.GetNumAtoms() == 0 else mol)
write_smiles(valid, 'output.smi')
open('failed_records.txt', 'w').write('\n'.join(map(str, failed)) + '\n')
```

The helper is `scripts/molecule_io.py`. SDF properties remain on parsed molecules; its read helpers skip failures, so track source records explicitly. Use `references/reading.md` for safe SMILES parsing and property preservation.

Done when converted records and failures are recorded and stereo checked.
