# SMARTS matching

```python
import json
from rdkit import Chem
from scripts.substructure_search import get_matches
smarts = open('query.smarts').readline().strip()
smiles = [x.strip() for x in open('library.smi') if x.strip()]
results = []
for i, s in enumerate(smiles):
    mol = Chem.MolFromSmiles(s)
    results.append({'index': i, 'smiles': s, 'matches': get_matches(mol, smarts) if mol else None})
json.dump(results, open('matches.json', 'w'), indent=2)
```

The helper is `scripts/substructure_search.py`; None denotes invalid SMILES. Match indices follow query atom order. For combined include/exclude filters, recursive SMARTS and stereochemistry see `references/smarts.md`, `references/matching.md` and `references/errors.md`.

Done when matches, input indices and invalid records are saved.
