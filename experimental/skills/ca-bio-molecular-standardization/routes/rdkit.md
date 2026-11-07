# Explicit RDKit standardization

```python
import json
from scripts.standardize_library import rdkit_standardize
smiles = [line.strip() for line in open('input.smi') if line.strip()]
json.dump([{'input': s, 'result': rdkit_standardize(s, keep_isotopes=False)} for s in smiles], open('standardized.json', 'w'), indent=2)
```

The default strips isotope labels, selects the largest fragment and canonicalizes tautomer; choose per-input isotope handling when needed. Tautomer choice is a representation, not a predicted dominant state. See `references/rdkit.md` and `references/representation.md`.

Done when outputs and charge, stereo, isotope and tautomer decisions are recorded.
