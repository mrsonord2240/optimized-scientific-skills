# Maximum common substructure

```python
import json
from rdkit import Chem
from scripts.similarity_search import find_mcs
smiles = [x.strip() for x in open('library.smi') if x.strip()]
molecules = [Chem.MolFromSmiles(s) for s in smiles]
if not molecules or any(m is None for m in molecules): raise ValueError('Invalid or empty input')
json.dump(find_mcs(molecules, timeout=60), open('mcs.json', 'w'), indent=2)
```

Inspect canceled before interpreting a small MCS; see `references/mcs-cliffs.md`.

Done when SMARTS, size and timeout status are recorded.
