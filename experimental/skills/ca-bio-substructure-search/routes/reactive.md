# Reactive-group flags

```python
import json
from rdkit import Chem
from scripts.reactive_alerts import reactive_filter
smiles = [x.strip() for x in open('library.smi') if x.strip()]
rows = []
for i, s in enumerate(smiles):
    mol = Chem.MolFromSmiles(s)
    rows.append({'index': i, 'smiles': s, 'alert': reactive_filter(mol) if mol else None})
json.dump(rows, open('reactive_flags.json', 'w'), indent=2)
```

The helper is `scripts/reactive_alerts.py`; it returns the first matching class, and None denotes invalid input. Treat broad alerts as flags in context; intended covalent warheads may be relevant chemistry. See `references/alerts.md`.

Done when flags and invalid inputs are saved for review.
