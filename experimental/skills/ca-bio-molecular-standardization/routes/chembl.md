# ChEMBL parent standardization

```python
import json
from scripts.standardize_library import chembl_standardize
smiles = [line.strip() for line in open('input.smi') if line.strip()]
json.dump([{'input': s, 'result': chembl_standardize(s)} for s in smiles], open('standardized.json', 'w'), indent=2)
```

The helper is `scripts/standardize_library.py`. Its default flags ambiguous multifragment parents; use the largest-fragment fallback only under an explicit policy. ChEMBL does not canonicalize tautomers here. See `references/chembl.md`, `references/policies.md` and `references/errors.md`.

Done when inputs, outputs and statuses are retained with the chosen policy.
